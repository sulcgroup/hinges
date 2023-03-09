# Leaf spring engine simulation analysis

## Raw data

The full simulation trajectories can be found in [this Google Drive](https://drive.google.com/drive/folders/1KoocIZRPcRJ7us0q695_ya3cXNxYUn2t)

## Extracted data

The simulation trajectories were analyzed using the get_k.py script found in [scripts](https://github.com/sulcgroup/hinges/tree/main/scripts) to extract the angles between the two arms of the engine as well as the distance between the two ends of the bridge.  The two datasets were essentially the same so angle between arms was used for all futher analysis as it is less susceptible to local fluctuations.  The angle data sets were extracted from the output of get_k.py using extract_data.py (also in [scripts](https://github.com/sulcgroup/hinges/tree/main/scripts)).  The processed angle data can be found in [data](https://github.com/sulcgroup/hinges/tree/main/data)

## Analysis

Most analysis including production of all figures found in the manuscript can be found in the notebook, final_figures.ipynb, found in [scripts](https://github.com/sulcgroup/hinges/tree/main/scripts)

Analysis of hydrogen bond patterns in the theoretically single-stranded flexure regions was done using hinge_bonds.py (also in [scripts](https://github.com/sulcgroup/hinges/tree/main/scripts)).

## Output

Output from the final_figure.ipynb script can be found in [final_plots](https://github.com/sulcgroup/hinges/tree/main/final_plots).
