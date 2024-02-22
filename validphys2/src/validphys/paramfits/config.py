"""
Configuration class for the paramfits module
"""
from collections.abc import Mapping, Sequence
import re

from reportengine.configparser import Config, ConfigError, element_of
from validphys import tableloader, utils
from validphys.loader import LoaderError


class ParamfitsConfig(Config):
    def _get_table(self, loader_func, fname, config_rel_path):
        # TODO: This is here because I am extremely unconvinced it is the
        # right interface. There should be something more specific at the
        # reportengine level. There is some undeniable ugliness in referencing
        # self.loader, which does not exist, but it is still better than
        # creating a base class in an isolated file to avoid the circular import.
        try:
            res = self.loader.check_vp_output_file(
                fname.strip(), extra_paths=['.', config_rel_path]
            )
        except LoaderError as e:
            raise ConfigError(e) from e

        try:
            df = loader_func(res)
        except Exception as e:
            raise ConfigError(e) from e
        return df

    # TODO: Get rid of this
    def produce_fits_pdf_config(self, fits):
        """DO NOT USE. For internal use only,"""
        return [self.produce_fitpdf(fit)['pdf'] for fit in fits]

    # TODO: Try to remove the loop from here
    def produce_fits_name(self, fits):
        """NOTE: EXPERIMENTAL.
        Return a list with the ids of the fits"""
        return [fit.name for fit in fits]

    # TODO: Try to remove the loop from here
    def produce_fits_as(self, fits_pdf_config):
        """NOTE: EXPERIMENTAL. Return the as value of the fits, reading
        it from the installed pdf"""
        return [pdf.alphas_mz for pdf in fits_pdf_config]

    # TODO: Try to remove the loop from here.
    def produce_fits_as_from_fitdeclarations(self, fitdeclarations):
        """NOTE: EXPERIMENTAL. A hack to obtain fits_as from the
        fitdeclarations, without having to
        download and inspect the actual fits."""
        alpha_pattern = r'NNPDF\d\d(?:_[a-z]+)*_as_(\d\d\d\d).*'
        res = []
        for fit in fitdeclarations:
            m = re.match(alpha_pattern, fit)
            if not m:
                raise ConfigError(
                    f"Couldn't match fit name {fit} to the " f"pattern {alpha_pattern!r}"
                )
            res.append(float(m.group(1)) / 1000)
        return {'fits_as': res}

    def produce_fits_name_from_fitdeclarations(self, fitdeclarations):
        """Inject the names from the ``fitdeclarations`` as the fit_names
        property"""
        # Cast nslist away
        return {'fits_name': list(fitdeclarations)}

    def parse_blacklist_datasets(self, datasets: list):
        return datasets

    def produce_combine_dataspecs_pseudoreplicas_as(
        self, dataspecs, how='min', blacklist_datasets=()
    ):
        if not isinstance(dataspecs, Sequence):
            raise ConfigError(
                "dataspecs should be a sequence of mappings, not " f"{type(dataspecs).__name__}"
            )
        if how != 'min':
            raise ConfigError("Only min is implemented at the moment")

        dfs = []
        fitnames = []
        for spec in dataspecs:
            if not isinstance(spec, Mapping):
                raise ConfigError(
                    "dataspecs should be a sequence of mappings, "
                    f" but {spec} is {type(spec).__name__}"
                )
            with self.set_context(ns=self._curr_ns.new_child(spec)):
                _, df = self.parse_from_(None, 'fits_computed_pseudoreplicas_chi2', write=False)
                _, asval = self.parse_from_(None, 'fits_as', write=False)
                _, namelist = self.parse_from_(None, 'fits_name', write=False)
                if not dfs:
                    firstas = asval
                elif asval != firstas:
                    raise ConfigError("Expecting all as values to be the same")
                dfs.append(df)
                fitnames.append(namelist)
        finalnames = [utils.common_prefix(*ns) + '__combined' for ns in zip(*fitnames)]
        res = tableloader.combine_pseudoreplica_tables(
            dfs, finalnames, blacklist_datasets=blacklist_datasets
        )

        return {'fits_computed_pseudoreplicas_chi2': res}

    # TODO: autogenerate functions like this
    def parse_experiments_covmat_output(self, fname: str, config_rel_path):
        """NOTE: THIS INTERFACE IS EXPERIMENTAL AND MIGHT CHANGE IN THE FUTURE.
        Process the output CSV table of the experiments_covmat action
        and return an equivalent dataframe"""
        df = self._get_table(tableloader.load_experiments_covmat, fname, config_rel_path)
        return {'experiments_covmat': df}

    # TODO: Move these to their own module when that's supported by reportengine
    def produce_fits_matched_pseudoreplicas_chi2_output(self, pseudoreplicafile: str, fits_name):
        """DEPRECATED. DO NOT USE."""
        import numpy as np
        import pandas as pd

        try:
            df = pd.read_csv(pseudoreplicafile, sep='\t', index_col=[0, 1], header=[0, 1])
        except Exception as e:
            raise ConfigError(f"Failed to load the table: {e}") from e

        # Require that the fits are matched so we filer out some that are not
        # interesting or broken.
        try:
            df = df[fits_name]
        except Exception as e:
            raise ConfigError(
                "Mismatch between fits provided and fits " f"in the table {pseudoreplicafile}:\n{e}"
            ) from e
        ndataindexer = df.columns.get_locs([slice(None), 'ndata'])
        lentest = lambda x: len(np.unique(x.dropna())) <= 1
        samelens = df.iloc[:, ndataindexer].apply(lentest, axis=1).all()
        if not samelens:
            raise ConfigError("Incorrect data: Expected all experiments to have the same length.")
        chindexer = df.columns.get_locs([slice(None), 'central_chi2'])
        df = df.iloc[:, chindexer]
        df = df.swaplevel(0, 1)
        # Have it the way the existing functions like
        newcols = df.columns.set_levels([df.columns.levels[0], ['chi2']])
        df.columns = newcols
        return df

    def parse_fits_computed_pseudoreplicas_chi2_output(self, fname: str, config_rel_path):
        """Return a namespace (mapping) with the output of
        ``fits_computed_pseudoreplicas_chi2_table`` as read from the specified
        filename. Use a {@with@} block to pass it to the providers.
        The fit names must be provided explicitly."""
        return self._get_table(
            tableloader.load_fits_computed_pseudoreplicas_chi2, fname, config_rel_path
        )

    def produce_use_fits_computed_pseudoreplicas_chi2_output(
        self, fits_computed_pseudoreplicas_chi2_output, fits_name
    ):
        """Select the columns of the input file matching the fits."""
        df = fits_computed_pseudoreplicas_chi2_output
        try:
            df = df[fits_name]
        except Exception as e:
            raise ConfigError(f"Could not select the fit names from the table: {e}") from e

        return {'fits_computed_pseudoreplicas_chi2': df}

    def produce_use_fits_computed_psedorreplicas_chi2_output(
        self, fits_computed_psedorreplicas_chi2_output, fits_name
    ):
        """Select the columns of the input file matching the fits.

        Note: this is a copy of ``produce_use_fits_computed_pseudoreplicas_chi2_output``.
        It is here so that `fits_computed_pseudoreplicas_chi2` gets assigned whether
        `fits_computed_pseudoreplicas_chi2_output` or `fits_computed_psedorreplicas_chi2_output`
        is specified in the runcard. This is to ensure that old runcards still work.
        """
        df = fits_computed_psedorreplicas_chi2_output
        try:
            df = df[fits_name]
        except Exception as e:
            raise ConfigError(f"Could not select the fit names from the table: {e}") from e

        return {'fits_computed_pseudoreplicas_chi2': df}

    @element_of('extra_sums')
    def parse_extra_sum(self, s: dict):
        keys = {'dataset_item', 'components'}
        if s.keys() != keys:
            d1 = s.keys() - keys
            d2 = keys - s.keys
            if d1:
                raise ConfigError(f'Unable to parse extra_sum: unrecognized keys: {d1}')
            if d2:
                raise ConfigError(
                    f'Unable to parse extra_sum. The following keys are required: {d2}'
                )
            raise RuntimeError()
        return s

    def produce_fits_matched_pseudoreplicas_chi2_by_experiment_and_dataset(
        self, fits_computed_pseudoreplicas_chi2, prepend_total: bool = True, extra_sums=None
    ):
        """Take the table returned by
        ``fits_matched_pseudoreplicas_chi2_output`` and break it down
        by experiment. If `preprend_total` is True, the sum over experiments
        will be included.

        This provides a namespace list with `suptitle`, `ndata` and
        `fits_replica_data_correlated`.

        """

        def get_ndata(df):
            val = df.index.get_level_values(2).unique()
            if len(val) != 1:
                raise ConfigError(f"Found different number " f"of points in {df.name}")
            return val[0]

        df = fits_computed_pseudoreplicas_chi2

        if prepend_total:
            s = df.loc[(slice(None), 'Total'), :].groupby(level=3).sum(min_count=1)
            ndata = (
                df.loc[(slice(None), 'Total'), :].groupby(level=0).apply(get_ndata).sum(min_count=1)
            )
            total = [
                {
                    'experiment_label': 'Total',
                    'by_dataset': [
                        {'fits_replica_data_correlated': s, 'suptitle': 'Total', 'ndata': ndata}
                    ],
                }
            ]
        else:
            total = []

        expres = []
        for exp, expdf in df.groupby(level=0):
            d = {'experiment_label': exp}
            by_dataset = d['by_dataset'] = []
            for ds, dsdf in expdf.groupby(level=1):
                ndata = dsdf.groupby(level=0).apply(get_ndata).sum()
                dsdf.index = dsdf.index.droplevel([0, 1, 2])

                if ds == 'Total':
                    if exp != 'Total':
                        ds = f'{exp} Total'
                    by_dataset.insert(
                        0, {'fits_replica_data_correlated': dsdf, 'suptitle': ds, 'ndata': ndata}
                    )
                else:
                    by_dataset.append(
                        {'fits_replica_data_correlated': dsdf, 'suptitle': ds, 'ndata': ndata}
                    )

            expres.append(d)

        if extra_sums:
            dss = {d['suptitle'] for l in [*total, *expres] for d in l['by_dataset']}
            for es in extra_sums:
                label = es['dataset_item']
                components = es['components']
                diff = set(components) - dss
                if diff:
                    bad_item = next(iter(diff))
                    raise ConfigError(f"Unrecognized elements in extra_sum: {diff}", bad_item, dss)

                sliced = tableloader.get_extrasum_slice(df, components)
                s = sliced.groupby(level=3).sum(min_count=1)
                ndata = sliced.groupby(level=[0, 1]).apply(get_ndata).sum()

                total.append(
                    {
                        'experiment_label': label,
                        'by_dataset': [
                            {'fits_replica_data_correlated': s, 'suptitle': label, 'ndata': ndata}
                        ],
                    }
                )

        return [*total, *expres]

    def _breakup_by_dataset_item(self, l, dataset_items):
        if dataset_items is None:
            return [{**expdict, **dsdict} for expdict in l for dsdict in expdict['by_dataset']]

        positions = {ds: pos for ds, pos in zip(dataset_items, range(len(dataset_items)))}
        # NOTE: If you want duplicates for some reason, you'll need to rewrite
        # this algorithm.
        if len(positions) != len(dataset_items):
            raise ConfigError("'dataset_items' cannot have duplicates")

        res = {}

        for expdict in l:
            for dsdict in expdict['by_dataset']:
                dsname = dsdict['suptitle']
                if dsname in positions:
                    res[positions[dsname]] = {**expdict, **dsdict}
                    del positions[dsname]
        if positions:
            raise ConfigError(f"Unrecognized dataset_items: {list(positions)}")
        return [res[index] for index in range(len(dataset_items))]

    def produce_fits_matched_pseudoreplicas_chi2_by_dataset_item(
        self,
        fits_matched_pseudoreplicas_chi2_by_experiment_and_dataset,
        dataset_items: (list, type(None)) = None,
    ):
        """Reorder, filter and flatten the result of
        fits_matched_pseudoreplicas_chi2_by_experiment_and_dataset with the
        dataset_items list. If it's not provided, this is equivalent to:
        fits_matched_pseudoreplicas_chi2_by_experiment_and_dataset::by_dataset
        Otherwise, the dictionaries will be returned in the order they appear
        in dataset_items, if they appear.
        """
        l = fits_matched_pseudoreplicas_chi2_by_experiment_and_dataset
        return self._breakup_by_dataset_item(l, dataset_items)

    def produce_matched_pseudoreplicas_for_total(
        self, fits_matched_pseudoreplicas_chi2_by_experiment_and_dataset
    ):
        """Like ``fits_matched_pseudoreplicas_chi2_by_dataset_item``, but
        restriction the ``dataset_item`` selection to "Total" exclusively."""
        res = self.produce_fits_matched_pseudoreplicas_chi2_by_dataset_item(
            fits_matched_pseudoreplicas_chi2_by_experiment_and_dataset, ['Total']
        )
        return res

    def produce_fits_replica_data_correlated_for_total(self, matched_pseudoreplicas_for_total):
        """Extract `fits_replica_data_correlated` from
        `matched_pseudoreplicas_for_total`. This is a hack that cannot be
        done efficiently with collect because of
        https://github.com/NNPDF/reportengine/issues/8."""
        return [matched_pseudoreplicas_for_total[0]['fits_replica_data_correlated']]

    def parse_fits_chi2_paramfits_output(self, fname: str, config_rel_path):
        """Load the output of ``fits_chi2_table`` adapted to suit the
        ``paramfits`` module. The fit names must be provided explicitly."""
        return self._get_table(tableloader.load_adapted_fits_chi2_table, fname, config_rel_path)

    def produce_use_fits_chi2_paramfits_output(self, fits_chi2_paramfits_output, fits_name):
        ndatatable, chi2table = fits_chi2_paramfits_output
        try:
            chi2table = chi2table[fits_name]
        except Exception as e:
            raise ConfigError(f"Could not select the fit names from the table: {e}") from e
        return {'adapted_fits_chi2_table': chi2table, 'ndatatable': ndatatable}

    def produce_fits_central_chi2_by_experiment_and_dataset(
        self, adapted_fits_chi2_table, ndatatable, prepend_total=True, extra_sums=None
    ):
        """Take the table returned by
        ``fits_matched_pseudoreplicas_chi2_output`` and break it down
        by experiment. If `preprend_total` is True, the sum over experiments
        will be included.

        This provides a namespace list with `suptilte` and
        `fits_replica_data_correlated`."""

        df = adapted_fits_chi2_table

        if prepend_total:
            s = df.sort_index().loc[(slice(None), 'Total'), :].sum(min_count=1)
            total = [
                {
                    'experiment_label': 'Total',
                    'by_dataset': [
                        {
                            'fits_total_chi2': s,
                            'suptitle': 'Total',
                            'ndata': ndatatable.loc[(slice(None), 'Total')].sum(),
                        }
                    ],
                }
            ]
        else:
            total = []
        expres = []
        for exp, expdf in df.groupby(level='experiment'):
            d = {'experiment_label': exp}
            by_dataset = d['by_dataset'] = []
            for ds, dsdf in expdf.groupby(level=1):
                dsdf.index = dsdf.index.droplevel([0])
                ndata = ndatatable[(exp, ds)]
                if ds == 'Total':
                    ds = f'{exp} Total'
                    by_dataset.insert(0, {'fits_total_chi2': dsdf, 'suptitle': ds, 'ndata': ndata})
                else:
                    by_dataset.append({'fits_total_chi2': dsdf, 'suptitle': ds, 'ndata': ndata})

            expres.append(d)

        if extra_sums:
            dss = {d['suptitle'] for l in [*total, *expres] for d in l['by_dataset']}
            for es in extra_sums:
                label = es['dataset_item']
                components = es['components']
                diff = set(components) - dss
                if diff:
                    bad_item = next(iter(diff))
                    raise ConfigError(f"Unrecognised element in extra sum: {diff}", bad_item, dss)

                sliced = tableloader.get_extrasum_slice(df, components)
                s = sliced.sum()
                ndata = tableloader.get_extrasum_slice(ndatatable, components).sum()
                total.append(
                    {
                        'experiment_label': label,
                        'by_dataset': [{'fits_total_chi2': s, 'suptitle': label, 'ndata': ndata}],
                    }
                )

        return [*total, *expres]

    def produce_fits_central_chi2_by_dataset_item(
        self, fits_central_chi2_by_experiment_and_dataset, dataset_items: (list, type(None)) = None
    ):
        """Reorder, filter and flatten the result of
        fits_central_chi2_by_experiment_and_dataset with the
        dataset_items list. If it's not provided, this is equivalent to:
        fits_central_chi2_by_experiment_and_dataset::by_dataset
        Otherwise, the dictionaries will be returned in the order they appear
        in dataset_items, if they appear.
        """
        l = fits_central_chi2_by_experiment_and_dataset
        return self._breakup_by_dataset_item(l, dataset_items)

    def produce_fits_central_chi2_for_total(self, fits_central_chi2_by_experiment_and_dataset):
        res = self.produce_fits_central_chi2_by_dataset_item(
            fits_central_chi2_by_experiment_and_dataset, ['Total']
        )
        return res

    # Define aliases for functions with spelling mistakes in their names which have now been corrected
    # Do this so that old runcards still work
    produce_combine_dataspecs_pseudorreplicas_as = produce_combine_dataspecs_pseudoreplicas_as
    produce_fits_matched_pseudorreplicas_chi2_output = (
        produce_fits_matched_pseudoreplicas_chi2_output
    )
    parse_fits_computed_psedorreplicas_chi2_output = parse_fits_computed_pseudoreplicas_chi2_output
    produce_fits_matched_pseudorreplicas_chi2_by_experiment_and_dataset = (
        produce_fits_matched_pseudoreplicas_chi2_by_experiment_and_dataset
    )
    produce_fits_matched_pseudorreplicas_chi2_by_dataset_item = (
        produce_fits_matched_pseudoreplicas_chi2_by_dataset_item
    )
    produce_matched_pseudorreplcias_for_total = produce_matched_pseudoreplicas_for_total
