%#############################################################
%
% Fichier  : ClassThesis.cls
% Auteur   : Erwann Camberlin (inspiré de Adrien Rouviere)
%
% Définition de la classe ClassThesis pour l'import des packages 
% et la mise en page du manuscrit de thèse
%
%#############################################################

\NeedsTeXFormat{LaTeX2e}                                        % Nécessite LaTeX2e ou une version plus récente
\ProvidesClass{0_Sources/ClassThesis}[2026/01/13]               % Déclare la classe "ClassThesis"

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}      % Transfère toutes les options non gérées à la classe "book"
\ProcessOptions \relax                                          % Traite les options déclarées et applique celles passées à la classe
\LoadClass[11pt, a4paper]{book}                                 % Charge la classe "book" pour du papier A4 et applique une police à 11pt

% ?-----------------------------------------------------------
% ? SILENCE
% ?-----------------------------------------------------------

\RequirePackage{silence}                                        % Permet de filtrer ou supprimer certains avertissements LaTeX
    \WarningFilter{latex}{You have requested}
    \WarningFilter{latex}{Command \under}

% ?-----------------------------------------------------------
% ? BOITE À OUTILS
% ?-----------------------------------------------------------

\RequirePackage{morewrites}                                     % Permet de résoudre l’erreur “No room for a new \write”
\RequirePackage{keyval}                                         % Permet de gérer des options ou paramètres sous forme clé=valeur dans les commandes ou packages
\RequirePackage{templatetools}                                  % Fournit des outils pour créer des templates ou des classes LaTeX
\RequirePackage{xargs}                                          % Fournit des versions étendues de \newcommand et de fonctions similaires
\RequirePackage{etoolbox}                                       % Fournit des outils de programmation avancés (conditions, boucles, hooks…)

% ?-----------------------------------------------------------
% ? INFOS SUR LA THÈSE
% ?-----------------------------------------------------------

\makeatletter
    \DeclareRobustCommand\Titre[1]{\gdef\@Titre{#1}}
    \newcommand{\theTitre}{\@Titre}                             % Commande stockant le titre de la thèse
    \DeclareRobustCommand\Auteur[1]{\gdef\@Auteur{#1}}
    \newcommand{\theAuteur}{\@Auteur}                           % Commande stockant le nom de l'auteur de la thèse
    \DeclareRobustCommand\Soutenance[1]{\gdef\@Soutenance{#1}}
    \newcommand{\theSoutenance}{\@Soutenance}                   % Commande stockant la date de soutenance de la thèse
\makeatother

% ?-----------------------------------------------------------
% ? CONFIGURATION & GÉOMETRIE
% ?-----------------------------------------------------------

\RequirePackage{geometry}                                       % Configuration de la mise en page
    \geometry{inner=2.5cm,
                outer=2cm,
                top=2.5cm,
                bottom=2.5cm,
                headheight=13.6pt}

\RequirePackage[indent]{parskip}                                % Gère l'espacement entre les paragraphes et l'indentation ('indent' conserve l'indentation du premier paragraphe)
\RequirePackage{lscape}                                         % Permet de mettre certaines pages en orientation paysage
\RequirePackage{setspace}                                       % Permet de gérer l'interligne : simple, 1,5 ou double (\singlespacing, \onehalfspacing, \doublespacing)
\RequirePackage{emptypage}                                      % Supprime les numéros et en-têtes sur les pages vides générées automatiquement (ex : \cleardoublepage)
\RequirePackage{layouts}                                        % Permet d'afficher les éléments de la mise en page (marges, colonnes, texte…) pour vérification
\RequirePackage[section]{placeins}                              % Place automatiquement un FloatBarrier (empêche les figures/tableaux de passer dans la section suivante) à chaque nouvelle section
\RequirePackage{multicol}                                       % Fournit un environnement pour gérer des contenus en plusieurs colonnes

% ?-----------------------------------------------------------
% ? LANGUE
% ?-----------------------------------------------------------

\RequirePackage[french]{babel}                                  % Gestion de la typographie française
    \addto\captionsfrench{                                      % Redéfinition de certains termes en français
        \def\tablename{Tableau}
        \renewcommand*\listfigurename{Liste des figures}
        \renewcommand*\listtablename{Liste des tableaux}
        \renewcommand*\listalgorithmcfname{Liste des algorithmes}
        }

% ?-----------------------------------------------------------
% ? FONTS
% ?-----------------------------------------------------------

\RequirePackage[utf8]{inputenc}                                 % Permet de convertir divers encodages standards en "language interne LaTeX"
\RequirePackage[T1]{fontenc}                                    % Permet l'encodage de la police
\RequirePackage{fourier}                                        % Charge la police Fourier
\RequirePackage[babel=true,kerning=true]{microtype}             % Améliore la typographie
\RequirePackage{textcomp}                                       % Fournit des symboles supplémentaires compatibles avec T1 (ex : €, °, ©)
\RequirePackage{lipsum}                                         % Génère du texte factice "Lorem Ipsum"
\RequirePackage[normalem]{ulem}                                 % Permet le texte barré (\sout) et souligné (\uline)
\RequirePackage{lettrine}                                       % Permet de créer des lettrines au début des chapitres ou paragraphes
    \setlength{\DefaultNindent}{0pt}

% ?-----------------------------------------------------------
% ? COULEURS
% ?-----------------------------------------------------------

\RequirePackage[table,x11names,svgnames,dvipsnames]{xcolor}     % Permet de définir des couleurs

% * Couleurs pdf
\definecolor{pdfurlcolor}{rgb}{0,0,1}
\definecolor{pdffilecolor}{rgb}{0.7,0,0}
\definecolor{pdflinkcolor}{rgb}{0,0,0.6}
\definecolor{pdfcitecolor}{rgb}{0.5,0.08,0.08}

% * Couleurs Julia
\definecolor{julia_red}{rgb}{0.796,0.235,0.2}       % 203 60 51
\definecolor{julia_blue}{rgb}{0.251,0.388,0.847}    % 64, 99, 216
\definecolor{julia_green}{rgb}{0.22,0.596,0.149}    % 56, 152, 38
\definecolor{julia_purple}{rgb}{0.584,0.345,0.698}  % 149, 98, 178

\definecolor{abscolor}{rgb}{0.05098, 0.31765, 0.50588}          %TODO

% ?-----------------------------------------------------------
% ? EN-TÊTE
% ?-----------------------------------------------------------

\RequirePackage{fancyhdr}                                       % Permet de personnaliser les en-têtes et pieds de page
    \renewcommand{\chaptermark}[1]{\markboth{\chaptername\ \thechapter.\ #1}{}}     % En-têtes : "Chapitre N. Titre"
    \renewcommand{\sectionmark}[1]{\markright{\thesection. #1}}                     % En-têtes : "Section N.N. Titre"
    \pagestyle{fancy}                                           % Activation du style "fancy"
    \fancyhf{}                                                  % Réinitialisation de tous les en-têtes et pieds de page
        \fancyhead[LE,RO]{\thepage}                             % Affiche le numéro de page à gauche sur les pages paires (Left Even) et à droite sur les pages impaires (Right Odd)
        \fancyhead[LO]{\nouppercase{\rightmark}}                % Affiche le titre de section à gauche sur les pages impaires
        \fancyhead[RE]{\nouppercase{\leftmark}}                 % Affiche le titre de chapitre à droite sur les pages paires
        %\renewcommand{\headrulewidth}{0pt}                     %TODO

% ?-----------------------------------------------------------
% ? NAVIGATION
% ?-----------------------------------------------------------

\usepackage[hyphens]{url}                                       % Permet d’inclure des URLs dans le texte

\RequirePackage[french,pagebackref=true,hyperfootnotes=true]{hyperref}              % Active les hyperliens dans le document PDF
    \hypersetup{linktocpage=true,                               % Seuls les numéros de page dans la TOC sont cliquables
                    colorlinks=true,                            % true : liens colorés; false : encadrés
                    linkcolor=pdflinkcolor,                     % Couleur des liens internes
                    citecolor=pdfcitecolor,                     % Couleur des liens vers la bibliographie
                    filecolor=pdffilecolor,                     % Couleur des liens vers des fichiers
                    urlcolor=pdfurlcolor,                       % Couleur des liens externes (URLs)
                    hypertexnames=true,                         % Noms uniques pour les ancres hyperref
                    pdfdisplaydoctitle=true,                    % Titre du PDF affiché dans le lecteur
                    bookmarksopen = true,                       % Ouvre tous les bookmarks dans le PDF
                    bookmarksopenlevel = 0,                     % Niveau jusqu’auquel les bookmarks sont ouverts
                    bookmarksnumbered = true,                   % Numérote les bookmarks
                    pdfnewwindow = true,                        % Ouvre les liens externes dans une nouvelle fenêtre
                    pdfpagelayout = OneColumn,                  % Affichage en une seule colonne (une seule page)
                    pdfstartview={XYZ null null 1},             % Vue initiale du PDF
                    pdfremotestartview={XYZ null null 1},
                    pdfview={XYZ null null 1}}

\Urlmuskip=0mu plus 1mu                                         % Ajuste l’espacement pour le retour à la ligne des URLs

\RequirePackage{bookmark}                                       % Optimise et organise les bookmarks pour le PDF, améliore hyperref

% ?-----------------------------------------------------------
% ? TABLE DES MATIERES
% ?-----------------------------------------------------------

\RequirePackage[titles]{tocloft}                                % Permet de personnaliser la ToC, ToF et ToT
\RequirePackage[french,nohints]{minitoc}                        % Permet de générer des mini-ToC au début de chaque chapitre
    \setcounter{minitocdepth}{2}                                % Profondeur du mini-ToC (chapter=0, section=1, subsection=2)
    \dominitoc                                                  % Active la génèration des mini-ToC

% ?-----------------------------------------------------------
% ? PERSONALISATION
% ?-----------------------------------------------------------

\RequirePackage{sectsty}                                        % Permet de personnaliser l'apparence des titres de sections, chapitres, etc.
    \allsectionsfont{\usefont{OT1}{phv}{b}{n}}                  % Définit la police de tous les titres de section : Helvetica (phv), gras (b), normal (n)
\RequirePackage[explicit]{titlesec}                             % Permet une personnalisation avancée des titres

\RequirePackage{epigraph}                                       % Permet d'insérer des citations/épigraphes au début des parties ou chapitres
    \setlength\epigraphwidth{.7\textwidth}

\RequirePackage[shortlabels]{enumitem}                          % Permet de personnaliser les listes

% ?-----------------------------------------------------------
% ? CHAPITRES
% ?-----------------------------------------------------------

\RequirePackage[Lenny]{fncychap}                                % Permet de personnaliser l’apparence des chapitres
    \ChNameVar{\fontsize{18}{20}\scshape\selectfont}            % Définit le style du mot "Chapitre"
    \ChNumVar{\normalfont\fontsize{60}{62}\selectfont}          % Définit le style du numéro du chapitre
    \ChTitleVar{\Huge\bfseries}                                 % Définit le style du titre du chapitre
    \makeatletter                                               % Modification des marges supérieurs des chapitres
        \patchcmd{\@makechapterhead}{\vspace*{50\p@}}{\vspace*{25\p@}}{}{}          % Réduit l’espace vertical avant les chapitres numérotés
        \patchcmd{\@makeschapterhead}{\vspace*{50\p@}}{\vspace*{25\p@}}{}{}         % Réduit l’espace vertical avant les chapitres non numérotés
    \makeatother

% ?-----------------------------------------------------------
% ? MATHEMATIQUES
% ?-----------------------------------------------------------

\RequirePackage{amsmath,mathtools}                              % Extensions AMS pour les mathématiques
    \allowdisplaybreaks                                         % Autorise les équations alignées à se couper sur plusieurs pages si nécessaire
\RequirePackage{amssymb}                                        % Fournit des symboles mathématiques supplémentaires

\RequirePackage[sfdefault=fav,scaled=0.9]{isomath}              % Permet de configurer la typographie des mathématiques
\RequirePackage{empheq}                                         % Permet de mettre en forme les environnements mathématiques

\newcommand\hmmax{0}
\newcommand\bmmax{0}
\RequirePackage{bm}                                             % Permet de mettre en gras des symboles mathématiques

\RequirePackage[fleqn,amsstyle]{cases}                          % Améliore l’environnement \begin{cases}

\RequirePackage{esint}                                          % Fournit des symboles d’intégrales plus jolis
\RequirePackage{nicefrac}                                       % Permet de créer des fractions "inline"

\RequirePackage[amsmath,amsthm,thmmarks,framed,hyperref,thref]{ntheorem}            % Permet de créer et personnaliser des environnements type théorème, proposition, définition
\RequirePackage{framed}                                         % Fournit des environnements encadrés

\RequirePackage[load-configurations = abbreviations]{siunitx}   % Gestion des unités SI
    \sisetup{unit-mode = text,
            locale = FR, detect-all,
            zero-decimal-to-integer,
            list-final-separator = { et },
            list-pair-separator = { et },
            range-phrase = { à }}
    \DeclareSIPostPower\quatre{4}

% ?-----------------------------------------------------------
% ? TABLEAUX
% ?-----------------------------------------------------------

\RequirePackage{multirow}                                       % Permet de fusionner plusieurs lignes dans un tableau
\RequirePackage{booktabs}                                       % Permet de créer des tableaux plus complexes
\RequirePackage{longtable}                                      % Permet de créer des tableaux qui s’étendent sur plusieurs pages
\RequirePackage[breakable, most]{tcolorbox}                     % Permet de créer des encadrés colorés et personnalisables

% ?-----------------------------------------------------------
% ? IMAGES
% ?-----------------------------------------------------------

\RequirePackage{graphicx}                                       % Permet d’inclure des images et des fichiers PDF dans le document
\RequirePackage{pdfpages}                                       % Permet d’inclure des pages entières de fichiers PDF
\RequirePackage{svg}                                            % Permet d’inclure des images SVG exportées depuis Inkscape

\RequirePackage{wrapfig}                                        % Permet d’insérer des figures flottantes à côté du texte

\RequirePackage[textfont={it},labelfont={bf}]{caption}          % Personnalisation des légendes par défaut (texte en italique, label en gras)
\RequirePackage[labelformat=simple]{subcaption}                 % Permet d'ajouter des légendes aux sous-figures
    \renewcommand\thesubfigure{(\alph{subfigure})}              % Format des sous-figures : (a), (b), ...

\RequirePackage{tikz}                                           % Permet de créer des dessins vectoriels avec TikZ
    \usetikzlibrary{angles}
    \usetikzlibrary{arrows,arrows.meta}
    \usetikzlibrary{calc}
    \usetikzlibrary{decorations.markings}
    \usetikzlibrary{patterns}
    \usetikzlibrary{plotmarks}
    \usetikzlibrary{shadows}
    \usetikzlibrary{snakes}
    \usetikzlibrary{spy}
    \usetikzlibrary{quotes}
    \usetikzlibrary{external}

\RequirePackage{pgf}                                            % Permet de créer des dessins vectoriels avec pgf
    \pgfkeys{/pgf/number format/.cd,1000 sep={}}

\RequirePackage{pgfplots}                                       % Permet de tracer des courbes et graphiques avec pgfplots
    \pgfplotsset{compat=newest}
    \pgfplotsset{plot coordinates/math parser=false}
    \usepgfplotslibrary{units}
    \usepgfplotslibrary{fillbetween}
    \usepgfplotslibrary{groupplots}
    \usepgfplotslibrary{external}

\RequirePackage{listofitems}                                    % Permet de lire et manipuler des listes d’items

\RequirePackage[outline]{contour}                               % Permet de créer des contours autour du texte
    \contourlength{1.4pt}

% ?-----------------------------------------------------------
% ? CODE
% ?-----------------------------------------------------------

\RequirePackage{listings}                                       % Permet d’insérer du code source avec coloration syntaxique
\RequirePackage{0_Sources/Pythonlistings}                       % Charge des configurations personnalisées pour insérer du code Python
\RequirePackage{0_Sources/Julialistings}                        % Charge des configurations personnalisées pour insérer du code Julia
\RequirePackage[autoload=true,linenumbers=true]{0_Sources/jlcode/jlcode}
\RequirePackage[numbered,framed]{matlab-prettifier}             % Permet d’insérer du code Matlab
\RequirePackage{verbatim}                                       % Permet d’insérer du texte brut ou du code TeX

\RequirePackage[ruled,vlined,linesnumbered,french,onelanguage]{algorithm2e}         % Permet de rédiger des algorithmes

% ?-----------------------------------------------------------
% ? NOMENCLATURE
% ?-----------------------------------------------------------

%\RequirePackage[automake,xindy={language=french},acronym,indexonlyfirst,nonumberlist]{glossaries-extra}        % Glossaire
\RequirePackage[automake,makeindex,acronym,indexonlyfirst,nonumberlist]{glossaries-extra}       % Permet de gérer la nomenclature

% ?-----------------------------------------------------------
% ? BIBLIOGRAPHIE
% ?-----------------------------------------------------------

\RequirePackage[sort&compress,comma,numbers]{natbib}            % Gestion de la bibliographie : citation numérique, citations multiples séparées par virgules, trie et compresse les citations consécutives
\RequirePackage{bibentry}
    \nobibliography*

%#############################################################
\endinput
%#############################################################