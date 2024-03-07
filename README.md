### How to build and execute program:
```
g++ -o main main.cpp
./main
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
