### Python Scripts for the research article "Deficiency of ZC3HC1 modulates vascular smooth muscle cell phenotype and increases neointima formation"

~ STRING DB network:
* xx_ZC3_string_db Python scripts creates STRING DB networks using the REST API and visualizes this network using networkx.

~ handling confocal images
* The script plot_confocal_IF.py can handle LEICA / .lif-files and convert it to png-files and creates super-projections from z-stacks
* cut_all_channels.py rotates and cuts the image at a defined region of interest (define parameters = shift, top, bottom)
* The script show_img_as_matrix_NIPA.py orders images as matrix (rows = channels/ staining) 
