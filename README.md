Intro
=====

Crée à partir des scripts de Pierre. L'application prend en compte 4 fichiers
en entrée. Elle applique une normalisation log(cpm) sur les données puis
crée diverses *Distributions* pour décrire celle-ci. (Onglet Descriptive
analysis) Peut être est-il possible d'incorporer par la suite la partie
écrite par Pierre dans un onfle Statistical analysis... ?

Input specifications
====================

-   Liste des gènes essentiels contenues dans une colonne au format
    texte, sans header.
-   Liste des gènes non-essentiels contenues dans une colonne au format
    texte, sans header.
-   Table de comptage au format csv. Les rownames correspondent aux noms
    de guides, les colnames aux échantillons.
-   Fichier de description des échantillons contentant une seule colonne
    au format texte. Respectant le format suivant.
    SAMPLENAME|Réplica-celltype-timepoint...

Contact
=======

-   <pierre.gestraud@curie.fr>
-   <clement.benoit@curie.fr>
