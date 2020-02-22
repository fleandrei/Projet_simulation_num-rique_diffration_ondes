

I) Ce projet contient:

   1) Dossier "code":
    Dossier contenant tous les fichier code du projet. Tous les fichiers présent    dans ce dossier ne sont pas tous utilisé dans l'exécutable final.
    Les fichiers utilisés par le programme sont:

	- diffraction_para.jl: Contient toutes les fonctions liées à la résolution du problème en lui même: construction de A (en séquentiel et parallèle), B, C, Dm, calcule de la valeur de l'onde en
	  un point précis...
	- polar.jl: Contient les fonctions relatives à la manipulation des coordonnées polaires.
	- initBoules.jl: Contient des fonctions pour générer des disques sur le plan

   2) Dossier "result":
    Dossier contenant des images obtenus par notre programme. Lorsqu'on lance le programme, c'est à cet endroit qu'est crée l'image. Le nom de l'image indique le nombre de disques ainsi que le mode
    d'exécution (séquentiel ou parallèle)

   3) Dossier global du projet:
     Fichier exécutable "run.jl" 


II) Lancer le programme:

    L'exécutable s'appelle "run.jl" et se trouve dans le dossier global du projet.
    Lancez le en ligne de commande avec julia et laissez vous guider par les instructions qui apparaissent dans le terminal.
    (Assurez vous que "julia" est bien la commande lancant julia, sinon il faut utiliser la commande qui lance julia)
    
   julia run.jl
    
    Il vous sera demandé si vous voulez faire apparaitre les disques de manière
    aléatoire (le nombre, la taille et la position des disques est aléatoire), ou de manière ordonnée en grille. Vous aurez alors la possibilité de choisir le nombre de disques, qui auront alors tous
    la même taille et seront placés les uns à coté des autres.
    On vous demandera également si vous voulez utiliser le parallèlisme ou non. Si vous choisissez de l'utiliser, alors le programme utilisera tous les coeurs disponible sur le PC.

    Par ailleurs, durant l'exécution, divers informations apparaîtrons sur le terminal telles que notamment, les temps de calcul de B, A, C, l'image. Pour avoir une bonne estimation du temps
    d'exécution d'une fonction sur julia, il faut lancer 2 fois le calcul. En effet, la première fois il y a un temps supplémentaire qui apparaît et qui est lié à la compilation. Le veritable temps
    d'exéctution est donc donné la seconde fois. Pour ne pas faire durer le programme trop longtemps, on ne lance qu'une seule fois le calcul. Le temps affiché n'est donc pas exactement le temps que
    met réelement la fontion mais en est relativement bien approché.

    Enfin, à la fin, le programme calcule le RCS du dommaine est créer un graphique dans le dossier results
	  
