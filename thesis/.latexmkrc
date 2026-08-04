# Svi generisani fajlovi (aux, log, xdv, pdf, ...) idu u build/
$out_dir = 'build';

# XeLaTeX: xdv -> pdf preko xdvipdfmx
$pdf_mode = 5;

# BibTeX radi u istom output direktorijumu
$bibtex_use = 2;

# Cursor/VS Code GUI nema /Library/TeX/texbin u PATH-u
$ENV{'PATH'} = '/Library/TeX/texbin:' . ($ENV{'PATH'} // '');
