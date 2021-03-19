.. _namespaces:

Namespaces
==========

In this section we explain a little bit about namespaces. Namespaces are used frequently in
`validphys` runcards because these runcards are based on `reportengine`, for which namespaces
are a central concept. 

A namespace is a stack of python dictionaries which are indexed by a tuple called the 
`namespace specification (nsspec)`. A user gives some inputs in terms of `fuzzyspecs` and
these are translated into `nnspecs. This is mostly an advanced internal implementation detail,
but it is important in order to understand how several features work. Also, the abstraction leaks
into user-facing 
