---
title: Data theory comparison for {@title dataset_input@} (cuts={@use_cuts@})
author: Zahari Kassabov
keywords: cms2ddy
...

{@with PTO@}
#{@ptoname@}({@theoryid@})
Absolute
---------

{@ plot_fancy @}


Normalized
----------

{@ datanorm plot_fancy@}

$\chi^2$
----

###Tables

{@ with pdfs @}
#### {@pdf@}
{@dataset_chi2_table @}
{@ endwith @}

### Replica distributions

{@pdfs plot_chi2dist @}

{@ endwith @}
