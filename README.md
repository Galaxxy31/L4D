# L4D - LaTeX 4 Doctorants

> *(ou Left 4 Dead, au choix)*

---

## Présentation

Bienvenue sur ce dépôt. Si vous êtes en phase de rédaction de votre thèse, ce template a été conçu pour vous accompagner dans la production de votre manuscrit. Si vous avez du temps avant la rédaction, vous pouvez également consulter ce dépôt pour réaliser une installation LaTeX sur VSCode.

L'objectif de ce dépôt est simple : vous fournir une **structure LaTeX prête à l'emploi** qui vous évite de perdre du temps sur la configuration technique. Vous vous concentrez sur le contenu scientifique, le template s'occupe du reste.

**Ce que ce dépôt vous propose :**
- Un squelette de thèse complet et fonctionnel
- Une configuration LaTeX pour VSCode
- Des exemples pour les éléments standards (bibliographie, nomenclature, annexes)

---

## Table des matières

- [Présentation](#présentation)
- [1. Pourquoi ce template](#1-pourquoi-ce-template)
- [2. Avant de commencer](#2-avant-de-commencer)
- [3. Extensions VSCode](#3-extensions-vscode)
- [4. Mise en place de LaTeX sur VSCode](#4-mise-en-place-de-latex-sur-vscode)
- [5. Utilisation du template](#5-utilisation-du-template)
- [6. Contribution](#6-contribution)
- [Ressources](#ressources)

---

## 1. Pourquoi ce template

### Le contexte

La rédaction d'une thèse représente un défi important : transformer plusieurs années de recherche en un manuscrit de **200 à 500 pages** qui doit être à la fois scientifiquement rigoureux et formellement irréprochable.

Les difficultés rencontrées sont souvent les mêmes :
- Structurer un document de cette taille
- Configurer LaTeX pour qu'il compile correctement
- Gérer une bibliographie complète
- Organiser les références croisées, la nomenclature, les annexes

Ce template répond à ces problématiques en fournissant une base solide et testée.

### Les avantages

| Besoin | Solution apportée |
|--------|-------------------|
| **Gain de temps** | La structure est déjà en place, vous pouvez commencer directement à écrire |
| **Qualité du rendu** | Mise en page propre conforme aux standards académiques |
| **Compilation fiable** | Configuration optimisée pour les longs documents (références, bibliographie, TikZ, etc.) |
| **Modularité** | Chaque chapitre est un fichier indépendant — modification sans impact sur le reste |
| **Réduction des erreurs** | Moins de configuration manuelle = moins de risques d'oubli |

### Approche modulaire

Le template utilise une **organisation modulaire** : chaque partie du manuscrit est contenue dans un fichier séparé. Cette approche présente plusieurs avantages :

1. **Compilation ciblée** : pendant la rédaction, compiler uniquement le chapitre en cours
2. **Travail parallèle** : possibilité de travailler sur différentes sections simultanément
3. **Maintenance simplifiée** : modifier un élément n'affecte pas l'ensemble du document

### Organisation du template

| Dossier | Rôle |
|---------|------|
| `0_Sources/` | Fichiers de configuration LaTeX (classes, styles) — à ne modifier que si vous savez ce que vous faites |
| `1_FrontMatter/` | Éléments préliminaires (couverture, remerciements, nomenclature) |
| `2_Chapitres/` | Corps du manuscrit |
| `3_BackMatter/` | Bibliographie, valorisation, résumé |
| `4_Annexes/` | Compléments techniques (démonstrations, algorithmes) |
| `5_Images/` | Figures, schémas, graphiques |
| `6_Datas/` | Données brutes, scripts de traitement |

---

## 2. Avant de commencer

### Prérequis

Pour utiliser ce template, vous aurez besoin de :

- **Un éditeur de texte** : VSCode est fortement recommandé. Le template de manuscrit a été conçu pour compiler et fonctionner correctement avec cet éditeur. Vous pouvez télécharger VSCode [ici](https://code.visualstudio.com/download). Notez que le template fonctionne également sur Overleaf si vous préférez cette solution.
- **MiKTeX à jour** : La DSI de l'ONERA a la fâcheuse tendance de ne pas tenir régulièrement à jour les logiciels. Dans la barre de recherche Windows, cherchez `MiKTeX Console`. Dans `Help → About MiKTeX console...`, vérifiez que votre version est égale ou supérieure à `4.12`. Si ce n'est pas le cas, vous devrez réaliser une demande de mise à jour auprès de la DSI.

### Recommandations

- Prenez le temps de lire ce README avant de commencer
- Faites une copie du template avant de modifier les fichiers de configuration
- Sauvegardez régulièrement votre travail

---

## 3. Extensions VSCode

### Obligatoires

| Extension | Description |
|-----------|-------------|
| **LaTeX** | Extension permettant à l'éditeur de reconnaitre les fichiers `.tex` comme des documents LaTeX et colorer les mots-clés. |
| **LaTeX Workshop** | Extension principale pour travailler avec LaTeX dans VSCode. Elle gère la compilation automatique du PDF, propose l'autocomplétion des commandes LaTeX, affiche un aperçu en temps réel, et signale les erreurs de compilation directement dans l'éditeur. Indispensable pour ce template. |
| **LaTeX Utilities** | Complète LaTeX Workshop avec des fonctionnalités supplémentaires comme la gestion des références croisées, la navigation entre les fichiers, et des outils de débogage avancés. |

### Recommandées

| Extension | Description |
|-----------|-------------|
| **Code Spell Checker** | Vérificateur d'orthographe pour détecter les fautes de frappe dans votre texte. |
| **French - Code Spell Checker** | Pack de langue français pour le vérificateur d'orthographe ci-dessus. |
| **Colorful Comments** | Ajoute des couleurs aux commentaires dans le code. Utilisée dans le template pour faire ressortir les différentes parties. |
| **Error Lens** | Affiche les erreurs et avertissements directement sur la ligne concernée, sans avoir à ouvrir le panneau d'erreurs. |
| **vscode-pets** | Peut être au final l'extension la plus importante à installer pour la rédaction du manuscrit. Une petite distraction pendant les longues compilations ! Des animaux virtuels qui parcourent votre éditeur. |

---

## 4. Mise en place de LaTeX sur VSCode

### Configuration du projet

Le dépôt inclut un fichier `.vscode/settings.json` qui contient toute la configuration nécessaire pour travailler efficacement avec LaTeX :

1. Ouvrez VSCode dans le dossier du projet
2. Allez dans `Fichier → Préférences → Paramètres` (ou `Ctrl` + `,`)
3. Cliquez sur l'icône en haut à droite pour ouvrir `settings.json`
4. Copiez-collez le contenu du fichier `.vscode/settings.json` du dépôt

### Contenu du fichier settings.json

Le fichier de configuration contient plusieurs sections :

#### Terminal
Configure Git Bash comme terminal par défaut sur Windows pour une meilleure compatibilité avec les commandes LaTeX.

#### Vérificateur d'orthographe (CSpell)
- Activation du vérificateur d'orthographe
- Langues configurées : anglais et français
- Mode de suggestion : correction rapide (quickFix)
- Autorise les mots composés

#### LaTeX Workshop - Outils de compilation
Définit les commandes disponibles pour la compilation :
- **pdflatex** : Compilation principale avec options pour l'externalisation TikZ (`--shell-escape`), la synchronisation (`-synctex=1`), et l'affichage des erreurs (`-file-line-error`)
- **bibtex** : Génération de la bibliographie

#### LaTeX Workshop - Recettes de compilation
Quatre méthodes de compilation prédéfinies :
| Recette | Description |
|---------|-------------|
| `pdflatex` | Compilation simple (pour tests rapides) |
| `pdflatex * 2` | Deux compilations successives (pour les tables de matières) |
| `pdflatex → bibtex → pdflatex * 2` | Compilation complète avec bibliographie |
| `Full Compilation` | Compilation complète avec 3 passes LaTeX finales (pour glossaire et le double référencement de la bibliographie) |

#### Comportement automatique
- **AutoBuild** : Désactivé pour éviter les compilations intempestives
- **AutoClean** : Nettoyage automatique des fichiers intermédiaires après une compilation réussie
- **Nettoyage** : Liste complète des fichiers à supprimer (.aux, .log, .toc, .gls, etc.)

#### Affichage PDF
- **Visionneur** : Ouvrir dans un nouvel onglet VSCode
- **Zoom** : Ajusté à la largeur de la page

### Vérification de la configuration

Pour vérifier que tout est correctement configuré :
1. Ouvrez un fichier `.tex`
2. Vous devriez voir l'icône LaTeX dans la barre d'activités
3. Sous l'option `Build LaTeX Project`, vous devriez voir les 4 recettes de compilation disponibles dans le menu

---

## 5. Utilisation du template

### 5.1 Récupérer le dépôt

1. Cloner le dépôt
```bash
git clone https://github.com/Galaxxy31/L4D.git
cd L4D
```

2. Récupérer un zip du dépôt via l'onglet `Code` puis `Download ZIP`

### 5.2 Structure du manuscrit

La structure du template est une structure un peu classique d'un manuscrit mais libre à vous de la modifier comme vous le souhaitez.

#### Front Matter (éléments préliminaires)

| Fichier | Description |
|---------|-------------|
| `1_Couverture_de_these.pdf` | Page de couverture (générez-la via ADUM) |
| `2_Seconde_de_de_couverture.tex` | Page de titre interne |
| `3_Remerciements.tex` | Remerciements |
| `4_Nomenclature.tex` | Fichier contenant la liste des symboles et abréviations |

#### Main Matter (corps du texte)

Le manuscrit est organisé en **3 parties** :

**Partie 1 : État de l'art**
- `1_Contexte.tex` : contexte de la recherche
- `2_Theorie.tex` : fondements théoriques
- `3_Problematique.tex` : problématique

**Partie 2 : Développements**
- `4_Physique.tex` : aspects physiques et modélisation
- `5_Numerique.tex` : méthodes numériques
- `6_Experimentale.tex` : dispositifs expérimentaux

**Partie 3 : Résultats et discussions**
- `7_Resultats.tex` : présentation des résultats
- `8_Discussion.tex` : analyse et interprétation

#### Back Matter

| Fichier | Description |
|---------|-------------|
| `Valorisation.tex` | Publications et communications |
| `Bibliographie.bib` | Entrées bibliographiques (BibTeX) |
| `Resume.tex` | Résumé (français/anglais) |

#### Annexes

- `A_Demonstration.tex` : démonstrations mathématiques
- `B_Algorithmes.tex` : algorithmes et codes

### 5.3 Bibliographie

La bibliographie est gérée via un fichier BibTeX (`3_BackMatter/Bibliographie.bib`). Exemple d'entrée :

```bibtex
@article{exemple2024,
  author = {Nom, Prénom},
  title = {Titre de l'article},
  journal = {Nom du Journal},
  year = {2024},
  volume = {1},
  pages = {10-25}
}
```

Citation dans le texte : `\cite{exemple2024}`

Vous pouvez générer un fichier `.bib` via votre outil de gestion bibliographique (Zotero par exemple).

### 5.4 Nomenclature

La nomenclature utilise le package `glossaries`. Définissez les entrées dans `1_FrontMatter/Nomenclature.tex` comme suit :

```latex
\newglossaryentry{reynolds}{
    type=roman,
    name={Re},
    description={Nombre de Reynolds},
    sort={Re}
}
```

Utilisation dans le texte : `\gls{reynolds}`

### 5.5 Compilation

Pour générer le PDF final, rien de bien compliqué puisque nous avons configuré l'extension `LaTeX Workshop` dans le `settings.json`. Cliquez simplement sur la recette `Full Compilation`.

**Astuce** : Utilisez et commentez les `\includeonly{}` dans `manuscrit.tex` pour compiler uniquement certains chapitres pendant la rédaction.

**Attention** : Pour pouvoir utiliser les recettes, pensez à bien faire le focus sur un fichier `.tex` sinon la recette ne se lancera pas.

### 5.6 Externalisation TikZ

Pour les documents avec de nombreux diagrammes TikZ, activez l'externalisation (ligne 46 dans `manuscrit.tex`) :

```latex
\tikzexternalenable%
```

Cela compile les figures une seule fois et les réutilise, ce qui réduit considérablement le temps de compilation.

---

## 6. Contribution

Les contributions sont les bienvenues. Pour proposer des améliorations :

1. Ouvrez une issue pour décrire le problème ou la suggestion
2. Fork le dépôt
3. Créez une branche pour votre modification
4. Soumettez une Pull Request

---

## Ressources

- [Documentation LaTeX officielle](https://www.latex-project.org/help/documentation/)
- [Overleaf Learn](https://www.overleaf.com/learn)
- [TeX Stack Exchange](https://tex.stackexchange.com/)
- [LaTeX cheatsheet](https://quickref.me/latex)

---
