### How to build and execute program:
```
make clean
make
./sfml-app
```
### Generation of graph diagram:
- graph_dot_file.dot gets generated from running the main program.
- Run this code to then use the *.dot file to generate a *.png version of it, which will contain the diagram:
    - If graphviz is not already installed, run this in the terminal:
        ```
        sudo apt install graphviz
        ```
    - Main code, in Linux terminal:
        ```
        dot -Tpng graph_dot_file.dot -o graph_dot_file.png
        ```

### To create .svg file using graphviz from the terminal:
```
dot -Tsvg graph_dot_file.dot -o graph_dot_file.svg
```