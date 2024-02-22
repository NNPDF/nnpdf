"""
    StopWatch module for computing the time performance of n3fit
"""

import time


def get_time():
    """Returns the cputime and walltime
    Note: only relative times make sense"""
    cpu_time = time.process_time()
    wall_time = time.time()
    return cpu_time, wall_time


class StopWatch:
    """
    This class works as a stopwatch, upon initialization it will register
    the initialization time as `start` and times can be register by running
    the `.register_times(tag)` method.

    When the stopwatchn is stopped (with the `.stop()` method) it will generate
    two dictionaries with the relative times between every register time and the
    starting point.
    """

    start_key = "start"

    def __init__(self):
        self._cputimes = {}
        self._walltimes = {}
        self.reference_list = []
        self.register_times(self.start_key)

    def get_times(self, tag=None):
        """Return a tuple with the `tag` time of the watch
        defaults to the starting time

        Parameters
        ----------
            `tag`
                if none, defaults to start_key

        Returns
        -------
            (tag cpu time, tag wall time)
        """
        if tag is None:
            tag = self.start_key
        start_cpu = self._cputimes[tag]
        start_wall = self._walltimes[tag]
        return (start_cpu, start_wall)

    def register_times(self, tag):
        """Register an event named `tag`"""
        cputime, walltime = get_time()
        self._cputimes[tag] = cputime
        self._walltimes[tag] = walltime

    def stop(self):
        """Stops the stopwatch and create the output dictionary

        Returns
        -------
            - `dict_out`: a dictionary containing two dictionaries
                with all relatives cpu and walltimes
        """
        end_cpu, end_wall = get_time()
        start_cpu = self._cputimes[self.start_key]
        start_wall = self._walltimes[self.start_key]
        cpu_out = {"Total": end_cpu - start_cpu}
        wall_out = {"Total": end_wall - start_wall}
        # Compute all differences relative to start
        for key, item in self._cputimes.items():
            cpu_out[key] = item - start_cpu
            wall_out[key] = self._walltimes[key] - start_wall
        # Compute all differences relative to references
        for key, ref in self.reference_list:
            f_cpu = self._cputimes[key]
            f_wall = self._walltimes[key]
            # If a required reference is not in the dict
            # dont' fail, the stopwatch should not be dangerous
            i_cpu = self._cputimes.get(ref, f_cpu)
            i_wall = self._walltimes.get(ref, f_wall)
            key_out = f"{ref}_to_{key}"
            cpu_out[key_out] = f_cpu - i_cpu
            wall_out[key_out] = f_wall - i_wall
        dict_out = {"walltime": wall_out, "cputime": cpu_out}
        return dict_out

    def register_ref(self, tag, reference):
        """Register an event named tag and register a request
        to compute also the time difference between this event and `reference`
        """
        self.register_times(tag)
        self.reference_list.append((tag, reference))
