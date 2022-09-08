.. _server:

Servers
=======

The NNPDF collaboration employs a storage server that host various data
files, meant for both public and internal consumption. It hosts the
following URLs:

-  https://data.nnpdf.science: Hosts **public** NNPDF data such as PDF
   fits, releases etc.
-  https://vp.nnpdf.science: Hosts the `validphys <vp-index>`__
   report and displays an index of all of the reports.
-  https://wiki.nnpdf.science: Hosts the github wiki version.
-  https://packages.nnpdf.science/: Hosts the ``conda`` binary packages.
-  https://docs.nnpdf.science/: Hosts this documentation.

SSH is used to interact with the server, as described in
`Access <#access>`__ below.

The NNPDF server is a virtual machine (VM) maintained by the Centro
Calcolo at the physics department of the University of Milan. The
machine has 2 CPUs, 4GB of RAM, 1 TB of disk and it is running CentOS7.

The full disk is backed up every week by the Centro Calcolo. We perform
every Sunday a ``rsync`` from the ``/home/nnpdf`` folder to the
``nnpdf@lxplus`` account at CERN.


.. _server-access:

Access
------

User access
~~~~~~~~~~~

The access to the server is provided by
``ssh``/:ref:`vp-upload <upload>` with the following restrictions:

-  ``ssh`` access to ``root`` is forbidden.
-  There is a shared ``nnpdf`` user with low privileges. In order to
   login the user must send his public ssh key (usually in
   ``~/.ssh/id_rsa.pub``) to SC. The ``nnpdf`` is not allowed to login
   with password.

The ``nnpdf`` user shares a common ``/home/nnpdf`` folder where all
NNPDF material is stored. Public access to data is available for all
files in the ``/home/nnpdf/WEB`` folder. The ``validphys`` reports are
stored in ``/home/nnpdf/validphys-reports`` and the wiki in
``/home/nnpdf/WEB/wiki``.

Access for continuous deployment tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `conda packages <conda>`__ as well as the documentation are
automatically uploaded to the server by the Continous Integration
service (Travis), through an user called ``dummy`` which has further
reduction in privileges (it uses the `rssh
shell <https://linux.die.net/man/1/rssh>`__) and it is only allowed to
run the ``scp`` command. An accepted private key is stored securely in
the `Travis configuration <travis-variables>`__. The packages are
uploaded to ``/home/nnpdf/packages``.

HTTP access
~~~~~~~~~~~

Tools such as `conda <conda>`__ and `vp-get <download>`__ require access
to private URLs, which are password-protected, using HTTP basic_auth.
The access is granted by a ``/.netrc`` file containing the user and
password for the relevant servers. The ``/.netrc`` file is typically
generated at `installation <conda>`__ time. It should look similar to

::

   machine vp.nnpdf.science
       login nnpdf
       password <PASSWORD>

   machine packages.nnpdf.science
       login nnpdf
       password <PASSWORD>

The relevant passwords can be found
`here <https://www.wiki.ed.ac.uk/pages/viewpage.action?pageId=292165461>`__.


.. _web-scripts:

Web scripts
-----------

Validphys2 interacts with the NNPDF server by `downloading
resources <download>`__ and `uploading results <upload>`__.

The server scripts live in the validphys2 repository under the
``serverscripts`` folder.

The server side infrastructure that makes this possible currently aims
to be minimalistic, although it may need to be expanded to a more robust
web application in time. At the moment, only thing that is done is
maintaining some index files (currently for theories, fits, reports and
LHAPDF sets) which essentially list the files in a given directory. The
indexes are regenerated automatically when their correspondent folders
are modified. This is achieved by waiting for changes using the Linux
``inotify`` API and the
`asynwatch <https://github.com/Zaharid/asyncwatch>`__ module. These
scripts are often controlled by `cron jobs <#cron-jobs>`__.

The report index is used to display a webpage indexing the reports. It
retrieves extra information from a ``meta.yaml`` file in the top level
output directory, and (with lower priority) by parsing an ``index.html``
page contained in the report folder. Properties like title, author and
tags are retrieved from the HTML header of this file, and are expected
to be in the same format that Pandoc would have used to write them when
``meta.yaml`` is passed as a input. To produce it, the most convenient
way is setting the ``main`` flag of a report, as described in [Uploading
the result].

Additionally information from the mailing list is added to the index
page. Specifically we query the list for links to validphys reports and
add links to the emails next to the entries of the reports that are
mentioned. This is achieved with the ``index-email.py`` script. It needs
some authentication credentials to access the mailing list. The password
is stored in a file called ``EMAIL_BOT_PASSWORD``, which is not tracked
by git. The script outputs two files in the root folder,
``email_mentions.json`` which should be used by other applications (such
as the report indexer) and ``seen_emails_cache.pkl``, which is there to
avoid downloading emails that are already indexes. These files need to
be deleted when the format of the index is updated.

The report index uses the `DataTables <https://datatables.net/>`__ JS
library. It provides filtering and sorting capabilities to the indexes
tables. The source file is:

::

   serverscripts/validphys-reports/index.html

in the validphys2 directory. It should be updated from time to time to
highlight the most interesting reports at a given moment. This can be
done by for example displaying in a separate table at the beginning the
reports marked with some keyword (for example ‘nnpdf31’).

The Makefile inside will synchronize them with the server.

The report indexing script generates thumbnails in the
``WEB/thumbnails`` which are then associated to each report. This is
done by looking at the image files inside the ``figures`` folder of each
uploaded report (see the source of the script for more details). It is
expected that the server redirects the requests for
``vp.nnpdf.science/thumbnails`` to this folder.

Cron jobs
---------

The following cron jobs are registered for the ``nnpdf`` user:

-  every day at 4 AM run the ``index-email.py`` script.
-  at every reboot run ``index-reports.py``, ``index-fits.py``,
   ``index-hyperscan.py``, ``index-packahes-public.sh`` and
   ``index-packages-private.sh``, which monitor continuously the
   respective folders and create indexes that can be used by various
   applications. The first two are homegrown scripts (see `Web
   Scripts <#web-scripts>`__) and the later two use
   `conda-index <https://docs.conda.io/projects/conda-build/en/latest/resources/commands/conda-index.html>`__.

The following cron jobs are registered for the ``root`` user:

-  perform backup of ``/home/nnpdf`` in lxplus every Saturday at noon.
-  perform a certbot renew every Monday.
-  reboot every Sunday at 6am (in order to use new kernels).
-  perform system update every day.

Web server Configuration
------------------------

We are using ``nginx`` as a lightweight and simple web server engine.
The ``nginx`` initial configuration depends on the linux distribution in
use. Usually debian packages provide a ready-to-go version where the
``/etc/nginx/nginx.conf`` is already set to work with server blocks
(subdomains).

Other distributions like CentOS7 requires more gymnastics, here some
tricks:

-  make sure the ``/home/nnpdf`` folder can be accessed by the ``nginx``
   user
-  folders served by ``nginx`` must have permission 755
-  create 2 folders in ``/etc/nginx``: ``sites-available`` and
   ``sites-enabled``.
-  in the ``/etc/nginx/nginx.conf`` file indicate the new include path
   with ``include /etc/nginx/sites-enabled/*;`` and remove all location
   statements.
-  for each server block create a new file in
   ``/etc/nginx/sites-available`` and build a symbolic link in
   ``/etc/nginx/sites-enabled``.
-  remember to perform a ``sudo service nginx restart`` or
   ``sudo nginx -s reload`` to update the server block configuration.

Finally, here an example of ``nginx`` configuration for the
``vp.nnpdf.science`` server block without ssl encryption:

::

   server {
       listen  80;
       listen [::]:80;
       server_name vp.nnpdf.science;

       root /home/nnpdf/validphys-reports;
       location / {
         try_files $uri $uri/ =404;
           auth_basic "Restricted";
           auth_basic_user_file /home/nnpdf/validphys-reports/.htpasswd;
       }

       location /thumbnails {
           alias /home/nnpdf/thumbnails;
           try_files $uri $uri/ =404;
           auth_basic "Restricted";
         auth_basic_user_file /home/nnpdf/validphys-reports/.htpasswd;
       }
   }

Some URLs are password protected using the HTTP ``basic_auth``
mechanism. This is implemented by setting the corresponding
configuration in nginx, as shown above (specifically with the
``auth_basic`` and ``auth_basic_user_file`` keys). The ``.htpasswd``
files mentioned in the configuration are generated with the ``htpasswd``
tool.

DNS
~~~

The domain is hosted by `Namecheap <https://namecheap.com>`__, which
also manages the DNS entries. For each subdomain there is an ``A``
record always pointing to the same server IP, currently 159.149.47.24.
The subdomains are then handled as described in `Web
server <#web-server>`__. For example, a DNS query for
``packages.nnpdf.science`` returns

::

    $ dig packages.nnpdf.science

   ; <<>> DiG 9.11.3-1ubuntu1.7-Ubuntu <<>> packages.nnpdf.science
   ;; global options: +cmd
   ;; Got answer:
   ;; ->>HEADER<<- opcode: QUERY, status: NOERROR, id: 26766
   ;; flags: qr rd ra; QUERY: 1, ANSWER: 1, AUTHORITY: 0, ADDITIONAL: 1

   ;; OPT PSEUDOSECTION:
   ; EDNS: version: 0, flags:; udp: 65494
   ;; QUESTION SECTION:
   ;packages.nnpdf.science.        IN  A

   ;; ANSWER SECTION:
   packages.nnpdf.science. 1799    IN  A   159.149.47.24

   ;; Query time: 170 msec
   ;; SERVER: 127.0.0.53#53(127.0.0.53)
   ;; WHEN: Tue May 28 14:26:53 BST 2019
   ;; MSG SIZE  rcvd: 67

SSL encryption
~~~~~~~~~~~~~~

SSL encription is provided by `Let’s
Encrypt <https://letsencrypt.org>`__. The certificates are created using
the ``certbot`` program with the ``nginx`` module.

In order to create new ssl certificates, first prepare the ``nginx``
server block configuration file and then run the interactive command:

::

   sudo certbot --nginx -d <domain>

This will ask you several questions, including if you would like to
automatically update the ``nginx`` server block file. We fully recommend
this approach.

The certificate is automatically renewed by a `cron job <#cron-jobs>`__.
