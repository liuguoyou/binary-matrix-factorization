%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDFTEX commands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\pdfoutput=1
\pdfimageresolution=300
\pdfcompresslevel=5
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% packages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\usepackage{graphicx}
\usepackage{amssymb}
\usepackage{xcolor}
\usepackage{xspace}

\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amsthm}

\usepackage{url} % Formatting urls
%\usepackage[colorlinks=true]{hyperref} % To follows url links
%\usepackage{balance}
%\usepackage[lined,ruled,linesnumbered]{algorithm2e}
\usepackage{subfig}
%
\theoremstyle{definition}
\newtheorem{theorem}{Theorem}
\newtheorem{proposition}{Proposition}
\newtheorem{lemma}{Lemma}
\newtheorem{corollary}{Corollary}
\newtheorem{definition}{Definition}
\newtheorem{example}{Example}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\definecolor{NewTextFG}{rgb}{0.0,0.1,0.8}
\newcommand{\fix}{\marginpar{FIX}}
\newcommand{\pendcite}[1]{{\color{red}[#1]}}
\newcommand{\note}[1]{\colorbox{yellow}{\scriptsize NOTE: #1}} % a la Acrobat
\newcommand{\pending}[1]{{\color{red}\bf PENDING: #1}}
\newcommand{\missing}[1]{{\color{red}\bf MISSING: #1}}
\newcommand{\newtext}[1]{{\color{NewTextFG}{#1}}\xspace}
\newcommand{\patibulo}[1]{{\color{red}{#1}}\xspace}
\newenvironment{textnote}[1]{\colorbox{yellow}{\scriptsize #1$\gg$}}{\colorbox{yellow}{\scriptsize $\ll$}\xspace}
\newcommand{\best}[1]{\textbf{\color{blue}#1}\xspace}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conventions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcommand{\iter}[1]{^{(#1)}}
\newcommand{\BlackBox}{\rule{1.5ex}{1.5ex}}  % end of proof
\newenvironment{proofi}{\par\noindent{\bf Proof: }}{\hfill\BlackBox\\}
%\newcommand{\opt}[1]{\ensuremath{\hat{#1}}}
\def\opt{\ensuremath{^{*}}}
\providecommand{\argmin}{\mathop{\textup{argmin}}}
\def\gradient{\nabla}
\def\hessian{\nabla^2}
\def\reals{\ensuremath{\mathbb{R}}}
\def\naturals{\ensuremath{\mathbb{N}}}
\newcommand{\mxn}{m{\times}n}
\renewcommand{\vec}[1]{\ensuremath{\mathbf{\MakeLowercase{#1}}}}
\newcommand{\mat}[1]{\ensuremath{\mathbf{\MakeUppercase{#1}}}}
\newcommand{\rv}[1]{\ensuremath{\MakeUppercase{#1}}}
\newcommand{\dpdf}[1]{\mathrm{\ensuremath{\MakeUppercase{#1}}}}
\newcommand{\cpdf}[1]{\mathrm{\ensuremath{\MakeLowercase{#1}}}}
\newcommand{\st}{\ensuremath{\quad\mathrm{s.t.}\quad}}
\newcommand{\norm}[1]{\ensuremath{\left\|#1\right\|}}
\newcommand{\quant}[1]{\ensuremath{\left[#1\right]}}
\newcommand{\support}[1]{\mathrm{supp}(#1)}
\newcommand{\rankf}[1]{\mathrm{rank}(#1)}
\def\rank{\mathrm{rank}}
\newcommand{\fun}[1]{\mathrm{#1}}
\newcommand{\cost}[1]{\ensuremath{\ell_{#1}}\xspace}
\newcommand{\abs}[1]{\ensuremath{\left|#1\right|}}
\newcommand{\setdef}[1]{\ensuremath{\left\{#1\right\}}}
\newcommand{\vecdef}[1]{\ensuremath{\left(#1\right)}}
\newcommand{\setspan}{\ensuremath{\mathrm{span}}}
\newcommand{\svec}[1]{_{[#1]}}
\def\transp{^\intercal}
\newcommand{\refeq}[1]{(\ref{#1})}
\def\med{\ensuremath{\mathrm{med}}}
\newcommand{\col}[1]{_{#1}}
\newcommand{\row}[1]{^{#1}}
\def\Ident{\mat{I}}
\def\Risk{\ensuremath{\mathcal{R}}}
\def\Loss{\ensuremath{\ell}}
\def\Gaussian{\ensuremath{\mathcal{N}}}
\newcommand{\GaussianPDF}[2][\sigma]{\ensuremath{\frac{1}{\sqrt{2\pi#1^2}}e^{-\frac{#2^1}{2#1^2}} }}
\def\Exponential{\ensuremath{\mathrm{Exp}}}
\def\Bernoulli{\ensuremath{\mathrm{Ber}}}
\def\Laplacian{\ensuremath{\mathrm{Lap}}}
\newcommand{\LaplacianPDF}[2][\theta]{\ensuremath{\frac{1}{2#1}e^{-\frac{|#2|}{#1}} }}
\def\Indicator{\ensuremath{\mathbf{1}}}
\def\LG{\ensuremath{\mathrm{LG}}}
\def\MOEG{\ensuremath{\mathrm{MOEG}}}
\def\MOE{\ensuremath{\mathrm{MOE}}}
\newcommand{\spaze}[1]{\ensuremath{\mathcal{#1}}}
\def\sgn{\mathrm{sgn}}
\newcommand{\minimize}{\ensuremath{\mathrm{minimize}\quad}}
\newcommand{\maximize}{\ensuremath{\mathrm{maximize}\quad}}
\def\defeq{:=}
\def\assign{\leftarrow}
\def\nin{\in\!\!\!/}
\def\MOE{\mathrm{MOE}}
\def\MOEG{\mathrm{MOEG}}
\newcommand{\inva}[1]{\left(#1\right)^{-1}}
\def\inv{^{-1}}
\newcommand{\pinva}[1]{\left(#1\right)^{\dagger}}
\def\pinv{^{\dagger}}
\def\model{M}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYMBOL TABLE -- WHAT IS USED FOR WHAT -- COPY THIS INTO TARGET FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   | lower normal            | upper normal            | lower bold              | upper bold              |  mathcal 
% A |                         |                         |                         |                         |
% B |                         |                         |                         |                         |
% C |                         |                         |                         |                         |
% D |                         |                         |                         |                         |
% E |                         |                         |                         |                         |
% F |                         |                         |                         |                         |
% G |                         |                         |                         |                         |
% H |                         |                         |                         |                         |
% I |                         |                         |                         |                         |
% J |                         |                         |                         |                         |
% K |                         |                         |                         |                         |
% L |                         |                         |                         |                         |
% M |                         |                         |                         |                         |
% N |                         |                         |                         |                         |
% O |                         |                         |                         |                         |
% P |                         |                         |                         |                         |
% Q |                         |                         |                         |                         |
% R |                         |                         |                         |                         |
% S |                         |                         |                         |                         |
% T |                         |                         |                         |                         |
% U |                         |                         |                         |                         |
% V |                         |                         |                         |                         |
% W |                         |                         |                         |                         |
% X |                         |                         |                         |                         |
% Y |                         |                         |                         |                         |
% Z |                         |                         |                         |                         |
%
%         | lower normal            | upper normal            | lower bold              | upper bold              |  
% Alpha   |                         |                         |                         |                         |
% Beta    |                         |                         |                         |                         |
% Gamma   |                         |                         |                         |                         |
% Delta   |                         |                         |                         |                         |
% Epsilon |                         |                         |                         |                         |
% Zeta    |                         |                         |                         |                         |
% Eta     |                         |                         |                         |                         |
% Theta   |                         |                         |                         |                         |
% Iota    |                         |                         |                         |                         |
% Kappa   |                         |                         |                         |                         |
% Lambda  |                         |                         |                         |                         |
% Mu      |                         |                         |                         |                         |
% Nu      |                         |                         |                         |                         |
% Ksi     |                         |                         |                         |                         |
% Omicron |                         |                         |                         |                         |
% Pi      |                         |                         |                         |                         |
% Rho     |                         |                         |                         |                         |
% Sigma   |                         |                         |                         |                         |
% Tau     |                         |                         |                         |                         |
% Upsilon |                         |                         |                         |                         |
% Phi     |                         |                         |                         |                         |
% Chi     |                         |                         |                         |                         |
% Psi     |                         |                         |                         |                         |
% Omega   |                         |                         |                         |                         |
%
%
% dictionary
%
\def\dict{d}
\def\dictm{\mat{\dict}}
\def\dictv{\vec{\dict}}
\def\atom{\dictv}

\def\dpred{b}
\def\dpredm{\mat{\dpred}}
\def\dpredv{\vec{\dpred}}

\def\data{y}
\def\datam{\mat{\data}}
\def\datav{\vec{\data}}
\def\datarv{\rv{\data}}

\def\err{e}
\def\errm{\mat{\err}}
\def\errv{\vec{\err}}
\def\errrv{\rv{\err}}
%\def\errrv{\epsilon}

\def\noiserv{N}

\def\coef{a}
\def\coefm{\mat{\coef}}
\def\coefv{\vec{\coef}}
\def\coefrv{\rv{\coef}}
%\def\coefrv{\alpha}

\def\supp{z}
\def\suppm{\mat{\supp}}
\def\suppv{\vec{\supp}}
\def\supprv{\rv{\supp}}
%\def\supprv{\zeta}

\def\sign{s}
\def\signm{\mat{\sign}}
\def\signv{\vec{\sign}}
\def\signrv{\rv{\sign}}
%\def\signrv{\phi}

\def\corr{g}
\def\corrm{\mat{\corr}}
\def\corrv{\vec{\corr}}
\def\corrrv{\rv{\corr}}

\def\canon{\omega}

\def\val{v}
\def\valm{\mat{\val}}
\def\valv{\vec{\val}}
\def\valrv{\rv{\val}}
%\def\valrv{\nu}

\def\aux{u}
\def\auxm{\mat{\aux}}
\def\auxv{\vec{\aux}}
\def\auxrv{\rv{\aux}}


\def\ndims{m} 
\def\natoms{p}
\def\nsamples{n}
\def\nclasses{c}

\def\datavspace{\reals^{\ndims}}
\def\coefvspace{\reals^{\natoms}}
\def\suppvspace{\setdef{0,1}^{\natoms}}
\def\signvspace{\setdef{-1,1}^{\natoms}}
\def\dictvspace{\reals^{\ndims}}
\def\coefmspace{\reals^{\natoms{\times}\nsamples}}
\def\dictmspace{\reals^{\ndims{\times}\natoms}}
\def\datamspace{\reals^{\ndims{\times}\nsamples}}

\def\di{i} % dimension index
\def\si{j} % sample index
\def\ai{k} % atom index

\def\supset{\Gamma}
\def\supsize{\norm{\coefv}_0}
\def\maxnz{\gamma}
\def\subspace{\mathbb{V}\xspace}
\def\bits{\,\mathrm{bits}}
\def\Huber{\mathcal{H}}
\def\huberpar{h}

\def\modelclass{\mathcal{M}}
