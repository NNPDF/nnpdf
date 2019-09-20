# Accessing NNPDF computer resources

NNPDF uses several systems for communication and code development. They are
described here briefly, together with the means to get access.

[SSH keys](https://en.wikipedia.org/wiki/Secure_Shell) are used for [git](git)
and [server](server) access. You can read about how to set them up somewhere
over the internet; for example
[here](https://www.digitalocean.com/community/tutorials/how-to-set-up-ssh-keys--2).

## Git repositories

*See [Git & Github](git)*

The repositories needed to access and develop the code are under the NNPDF
Github organization

<https://github.com/NNPDF/nnpdf>

In oder to gain access, you need to set up a Github account and ask Stefano
Carrazza or Zahari Kassabov to include you in the list of members of the Github
organization.

## Mailing list

Emails are send to a [mailing
list](https://lists.cam.ac.uk/mailman/listinfo/ucam-nnpdf) provided by the
University of Cambridge. Users can subscribe to it filling the corresponding
form at the address


<https://lists.cam.ac.uk/mailman/listinfo/ucam-nnpdf>

You need to wait for some of the administrators (Maria Ubiali or Zahari
Kassabov) to approve your request.

Various mailing options can be changed from this URL

<https://lists.cam.ac.uk/mailman/options/ucam-nnpdf>

## NNPDF server

*See [Server Configuration](server)*

We maintain a storage server, which contains public and private data. Read
access to the private data (such as [validphys reports](vp-index)) is guarded by
an HTTP password. The relevant credentials are
[here](https://www.wiki.ed.ac.uk/pages/viewpage.action?pageId=292165461) (for
which you need access to the [wiki](#edinburgh-wiki), or else to ask someone).

Write access (such as to upload reports) is provided by SSH. You need to send
your public key to Stefano Carrazza or Zahari Kassabov.


## Edinburgh wiki

We uses wiki system hosted by the University of Edinburgh. It is accessible from

<https://www.wiki.ed.ac.uk/pages/viewpage.action?spaceKey=nnpdfwiki&title=The+NNPDF+Collaboration+wiki>

It stores data on historical archives and various organizational and strategy
documents (such as the list of [data to be
implemented](https://www.wiki.ed.ac.uk/display/nnpdfwiki/Experimental+data+and+applgrids+for+NNPDF4.0)),
In order to gain access, you need to ask Luigi del Debbio.

```eval_rst

.. warning::
    The Edinburgh wiki is not directly controlled by us, not backed up by us and
    it uses a closed system that makes it nearly impossible to export resources
    uploaded there should we be unable to use it in the future for any reason.
    Unlike emails and git repositories, no local copy of the arhive is easily
    available.  It is **HIGHLY DISCOURAGED** to upload there anything you care
    about losing at some point.  Instead, upload the relevant data to the NNPDF
    server and post links on the wiki if you are asked to do so.
```
