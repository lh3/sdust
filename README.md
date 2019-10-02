Sdust is a reimplementation of the [symmetric DUST algorithm][paper] for
finding low-complexity regions in DNA sequences. It gives identical output to
[NCBI's dustmasker][dm] except in assembly gaps, and is four times as fast. The
source code here was initially written for [minimap][mm] and later
[minimap2][mm2]. This repo is a standalone copy.

[paper]: http://www.ncbi.nlm.nih.gov/pubmed/16796549
[dm]: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker/
[mm]: https://github.com/lh3/minimap
[mm2]: https://github.com/lh3/minimap2
