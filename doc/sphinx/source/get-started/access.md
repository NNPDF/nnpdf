# Accessing NNPDF computer resources

NNPDF uses several systems for communication and code development. They are
described here briefly, together with the means to get access to them.

[SSH keys](https://en.wikipedia.org/wiki/Secure_Shell) are used for [git](git)
and [server](server) access. You can read about how to set them up somewhere
over the internet; for example
[here](https://www.digitalocean.com/community/tutorials/how-to-set-up-ssh-keys--2).

## Git repositories

*See [Git & GitHub](git)*

The repositories needed to access and develop the code are under the NNPDF
GitHub organization

<https://github.com/NNPDF/nnpdf>

In order to gain access, you need to set up a GitHub account and ask Stefano
Carrazza or Zahari Kassabov to include you in the list of members of the GitHub
organization.

```eval_rst
.. _mail:
```
## Mailing list

Emails are sent to a [mailing
list](https://lists.cam.ac.uk/mailman/listinfo/ucam-nnpdf) provided by the
University of Cambridge. Users can subscribe to it by filling in the corresponding
form at the address


<https://lists.cam.ac.uk/mailman/listinfo/ucam-nnpdf>

You need to wait for one of the administrators (Maria Ubiali or Zahari
Kassabov) to approve your request.

Various mailing options can be changed from this URL

<https://lists.cam.ac.uk/mailman/options/ucam-nnpdf>

## Rules for NNPDF email exchange

Never attach files, always send links.
Never reply including the original message.
Ensure that the message corresponds to a given specific thread:

* [NNPDF] For general NNPDF issues 
* [NNPDF_news] For matters related to NNPDF activities such as talks, conferences, meetings, grants, etc
* [NNPDF4_data] For matters related to the inclusion of new data in view of the NNPDF4.0 PDF set
* [N3FIT] For matters related to the new fitting code

Dedicated NNPDF threads for specific Physics projects, which might include non-NNPDF members

* [NNPDF_smallx] For small-x resummation: NNPDF + Marco Bonvini and Simone Marzani
* [NNFF] For the NNPDF fits of Fragmentation Functions
* [nNNPDF] For the NNPDF fits of nuclear parton distributions
* [NNPDF_th] For studies and fits related to the implementation of theory uncertainties
* [NNPDF_jets] For studies related to choice of jet observables and scales

## NNPDF server

*See [Server Configuration](server)*

We maintain a storage server, which contains public and private data. Read
access to the private data (such as [validphys reports](vp-index)) is guarded by
an HTTP password. The relevant credentials are
[here](https://www.wiki.ed.ac.uk/pages/viewpage.action?pageId=292165461) (for
which you need access to the [wiki](#edinburgh-wiki), or else to ask someone).

Write access (such as to [upload reports](upload)) is provided by SSH.
You need to send your public key to Stefano Carrazza or Zahari
Kassabov.


## Edinburgh wiki

We use a wiki system hosted by the University of Edinburgh. It is accessible from

<https://www.wiki.ed.ac.uk/pages/viewpage.action?spaceKey=nnpdfwiki&title=The+NNPDF+Collaboration+wiki>

It stores data on historical archives and various organizational and strategy
documents (such as the list of [data to be
implemented](https://www.wiki.ed.ac.uk/display/nnpdfwiki/Experimental+data+and+applgrids+for+NNPDF4.0)).
In order to gain access, you need to ask Luigi del Debbio.

```eval_rst

.. warning::
    The Edinburgh wiki is not directly controlled by us, not backed up by us and
    it uses a closed system that makes it nearly impossible to export resources
    uploaded there should we be unable to use it in the future for any reason.
    Unlike emails and git repositories, no local copy of the archive is easily
    available.  It is **HIGHLY DISCOURAGED** to upload anything to the wiki that you care
    about losing at some point.  Instead, :ref:`upload <upload>` the relevant data to the NNPDF
    server and post links on the wiki if you are asked to do so.
```
