# How to list the available resources

```{eval-rst}
.. _vp-list:
```

## Using `vp-list`

In order to check what resources are available locally and for download, use
`vp-list` which will print out the names of resources.

```bash
vp-list <resource type>
```

The options for resource type can be seen with `vp-list --help`.

```bash
$ vp-list --help
usage: vp-list [-h] [-r | -l] resource

vp-list Script which lists available resources locally and remotely

positional arguments:
  resource           The type of resource to check availability for (locally
                     and/or remotely). Choose from: theories, fits, pdfs,
                     datasets.

```

You can use the options `-l/--local-only` or `-r/--remote-only` to only check
for resources available locally or remotely respectively.

## Manually checking server - example with fits

You can also check manually on the storage servers for these resources. For example,
the bulk of the available fits can be found by going to the fits folder of the
NNPDF data server. Some other fits may be found in standalone folders in the home
folder of this server, but of course finding a specific fit here may require some
digging. For help in accessing the server,
please see [here](NNPDF-server). For information on how to download fits and
other resources,
please see the [Downloading resources](download) section of the vp-guide.