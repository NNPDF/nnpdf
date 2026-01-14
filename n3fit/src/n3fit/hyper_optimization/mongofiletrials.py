"""
Hyperopt trial object for parallel hyperoptimization with MongoDB.
Data are fetched from MongoDB databases and stored within the nnfit folder.

The workflow when running parallel hyperopt with mongodb is as follows:
    1. Submit the main "server" job. This can be submitted as part of n3fit or as a separate job,
    e.g., in the login node, running only the database.
    The main server job can run itself a worker job.
    2. Submit the worker jobs. The workers each run a single job.

Each of the workers need to run with the --parallel-hyperopt flag and --hyperopt <number of trials>,
where the number of trials is per-worker (so, for 1000 each, 2 workers, will give 2000 trials total).
Each job will first try to connect to the database and, failing that, will try to start the database.

The database is stored in the folder for the replica 1 (together with the trials) so that it gets
stored with vp-upload, allowing for a continuation of the fit.
"""

import json
import logging
import os
from pathlib import Path
import platform
import subprocess
import time

from bson import SON, ObjectId
from hyperopt.mongoexp import MongoTrials

from n3fit.backends import get_physical_gpus
from n3fit.hyper_optimization.filetrials import space_eval_trial

log = logging.getLogger(__name__)


def convert_bson_to_dict(obj):
    """
    Recursively convert a BSON object to a standard Python dictionary.

    This function is particularly useful for converting MongoDB query results,
    which may contain BSON types like ObjectId and SON, into a more manageable
    dictionary format.

    Parameters
    ----------
    obj : dict or bson.SON or list or any
        The object to convert. Can be a BSON object (like SON), a dictionary
        containing BSON types, a list of such objects, or any other type.

    Returns
    -------
    dict or list or any
        A Python dictionary with all BSON types converted to standard Python
        types (e.g., ObjectId converted to string). If the input is a list,
        returns a list of converted elements. For other types, returns the
        object as is.

    Examples
    --------
    >>> from bson import ObjectId, SON
    >>> sample_son = SON([('_id', ObjectId('507f1f77bcf86cd799439011')), ('name', 'John Doe')])
    >>> convert_bson_to_dict(sample_son)
    {'_id': '507f1f77bcf86cd799439011', 'name': 'John Doe'}

    >>> sample_list = [SON([('_id', ObjectId('507f1f77bcf86cd799439011')), ('name', 'John Doe')]), {'age': 30}]
    >>> convert_bson_to_dict(sample_list)
    [{'_id': '507f1f77bcf86cd799439011', 'name': 'John Doe'}, {'age': 30}]
    """
    if isinstance(obj, (SON, dict)):
        return {k: convert_bson_to_dict(v) for k, v in obj.items()}
    if isinstance(obj, ObjectId):
        return str(obj)  # or just return None if you don't need the ObjectId
    if isinstance(obj, list):
        return [convert_bson_to_dict(v) for v in obj]
    return obj


class MongodRunner:
    """Class to manage a MongoDB instance.

    This class is responsible for automatically creating and managing a MongoDB database
    using the `mongod` command. It allows for starting and stopping a MongoDB instance
    programmatically.

    Parameters
    ----------
        db_port: int
            MongoDB database connection port. Defaults to 27017.
        db_host: str
            hostname of the database
        db_name: str
            MongoDB database name. Defaults to "hyperopt-db".
    """

    def __init__(self, db_path="hyperopt-db", db_host=None, db_port=27017):
        self.db_path = db_path
        self.db_port = db_port
        self.db_host = db_host
        self._runner_job = None

    def is_up(self):
        """Checks whether the database is up."""
        from pymongo import MongoClient
        from pymongo.errors import ServerSelectionTimeoutError

        # If the db doesn't exist, and we didn't get a target to find it, no need to look further
        if not self.db_path.exists() and self.db_host is None:
            return False

        # If the db exists, check whether we know how to connect:
        if self.db_host is None:
            if (hostfile := self.db_path.with_suffix(".hostname")).exists():
                possible_host = hostfile.read_text()
        else:
            possible_host = self.db_host

        mc = MongoClient(host=possible_host, port=self.db_port, serverSelectionTimeoutMS=100)
        try:
            mc.admin.command("ismaster")
            self.db_host = possible_host
            log.info(f"Database up at {self.db_host}:{self.db_port}")
        except ServerSelectionTimeoutError:
            mc.close()
            return False
        return True

    def start(self):
        """Starts the MongoDB instance via `mongod` command."""
        args = [
            "mongod",
            "-quiet",
            "--wiredTigerCacheSizeGB",
            "16",
            "--dbpath",
            self.db_path,
            "--port",
            str(self.db_port),
            "--bind_ip_all",
        ]
        try:
            self.db_path.mkdir(exist_ok=True, parents=True)
            mongod = subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            cmd_str = " ".join([str(i) for i in args])
            # The job starting the database gets to write the hostname down
            # unless it has been set explicitly, in that case trust the user and use that
            if self.db_host is None:
                self.db_host = platform.node()
            self.db_path.with_suffix(".hostname").write_text(self.db_host)
            log.info(f"Started MongoDB database at {self.db_host}:{self.db_path} with {cmd_str}")
            self._runner_job = mongod
        except OSError as err:
            msg = f"Failed to execute {args}. Make sure you have MongoDB installed."
            self.db_path.with_suffix(".hostname").unlink(missing_ok=True)
            raise EnvironmentError(msg) from err

    def stop(self):
        """Stops `mongod` command."""
        if self._runner_job is None:
            return
        try:
            self._runner_job.terminate()
            self._runner_job.wait()
            log.info("Stopped mongod")
        except Exception as err:
            log.error(f"Failed to stop mongod: {err}")

    def __del__(self):
        # Change it to a context manager
        self.stop()


class MongoFileTrials(MongoTrials):
    """
    MongoDB implementation of :class:`n3fit.hyper_optimization.filetrials.FileTrials`.

    Parameters
    ----------
        replica_path: path
            Replica folder as generated by n3fit.
        mongod_runner: :py:class:MongodRunner
            Instance of MongodRunner with parameters such as db_host or db_port defining the database
        num_workers: int
            Number of MongoDB workers to be initiated concurrently. Defaults to 1.
        parameters: dict
            Dictionary of parameters on which we are doing hyperoptimization. Default to None.
        store_trial: bool
            If True, store data into json file. Default to True.
    """

    def __init__(
        self, replica_path, mongod_runner, *args, num_workers=1, parameters=None, **kwargs
    ):
        self.num_workers = num_workers

        # Define the connection string
        db_name = Path(mongod_runner.db_path).name
        host = mongod_runner.db_host
        port = mongod_runner.db_port
        self._mongo_str = f"{host}:{port}/{db_name}"

        self.workers = []
        self._runname = replica_path.parts[-3]  # corresponds to the runcard/runfolder name

        self._store_trial = False
        self._json_file = replica_path / "tries.json"
        self._parameters = parameters
        self._rstate = None
        self._dynamic_trials = []

        super().__init__(f"mongo://{self._mongo_str}/jobs", *args, **kwargs)

    @property
    def rstate(self):
        """Returns the rstate attribute; see :class:`n3fit.hyper_optimization.filetrials.FileTrials`."""
        return self._rstate

    @rstate.setter
    def rstate(self, random_generator):
        """Sets the rstate attribute; see :class:`n3fit.hyper_optimization.filetrials.FileTrials`."""
        self._rstate = random_generator

    def _set_dynamic_trials(self):
        """Converts self._trials to a dictionary and stores it in self._dynamic_trials."""
        self._dynamic_trials = [convert_bson_to_dict(item) for item in self._trials]

    def refresh(self):
        """Fetches data from mongo database and save to a json file."""
        super().refresh()

        # convert BSON object to a dictionary
        self._set_dynamic_trials()

        # write json to disk
        if self._store_trial:
            local_trials = []
            for idx, t in enumerate(self._dynamic_trials):
                local_trials.append(t)
                local_trials[idx]["misc"]["space_vals"] = space_eval_trial(self._parameters, t)

            all_to_str = json.dumps(local_trials, default=str)
            with open(self._json_file, "w") as f:
                f.write(all_to_str)

    # like in `FileTrials` the two methods below are implemented to avoid writing to the database twice
    def new_trial_ids(self, n):
        self._store_trial = False
        return super().new_trial_ids(n)

    def _insert_trial_docs(self, docs):
        self._store_trial = True
        return super()._insert_trial_docs(docs)

    def start_mongo_workers(
        self,
        workdir=None,
        exp_key=None,
        poll_interval=0.1,
        no_subprocesses=False,
        max_consecutive_failures=10,
        reserve_timeout=600,
    ):
        """Initiates all mongo workers simultaneously."""
        # get the number of gpu cards, if any
        gpus_all_physical_list = get_physical_gpus()
        num_gpus_available = len(gpus_all_physical_list)
        if not num_gpus_available:
            log.warning("No GPUs found in the system.")

        # construct the command to start a hyperopt-mongo-worker
        args = ["hyperopt-mongo-worker", "--mongo", self._mongo_str]
        if workdir:
            args.extend(["--workdir", workdir])
        if exp_key:
            args.extend(["--exp-key", exp_key])
        if poll_interval:
            args.extend(["--poll-interval", str(poll_interval)])
        if max_consecutive_failures:
            args.extend(["--max-consecutive-failures", str(max_consecutive_failures)])
        if reserve_timeout:
            args.extend(["--reserve-timeout", str(reserve_timeout)])
        if no_subprocesses:
            args.append("--no-subprocesses")

        # start the worker as a subprocess
        # TODO FIX FIXME
        # The mongotrials should not have to know anything about the outside world, each run will use the GPU it has access to and it is on the user to select the right one
        try:
            my_env = os.environ.copy()
            my_env["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"

            # create log files to redirect the mongo-workers output
            mongo_workers_logfile = f"mongo-worker_{self._runname}_{time.time()}.log"
            with open(mongo_workers_logfile, mode='w', encoding="utf-8") as log_file:
                # run mongo workers
                worker = subprocess.Popen(
                    args, env=my_env, stdout=log_file, stderr=subprocess.STDOUT
                )
                self.workers.append(worker)
            log.info("Started mongo worker")
        except OSError as err:
            msg = f"Failed to execute {args}. Make sure you have MongoDB installed."
            raise EnvironmentError(msg) from err

    def stop_mongo_workers(self):
        """Terminates all active mongo workers."""
        for worker in self.workers:
            try:
                worker.terminate()
                worker.wait()
                log.info(f"Stopped mongo worker {self.workers.index(worker)+1}/{self.num_workers}")
            except Exception as err:
                log.error(
                    f"Failed to stop mongo worker {self.workers.index(worker)+1}/{self.num_workers}: {err}"
                )
