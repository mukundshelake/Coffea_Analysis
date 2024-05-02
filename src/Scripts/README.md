flowchart TD
    %% ^ These subgraphs are identical, except for the links to them:

    %% Link *to* subgraph1: subgraph1 direction is maintained
    id1[(Data)] --> Skimmer --> id2{{skimmerOutput}} -->Extractor

    Extractor --> id3{{extractorOutput}}  -->Plotter
    %% Link *within* subgraph2:
    %% subgraph2 inherits the direction of the top-level graph (LR)
    unfilledTex --> Jinja --> filledReport

    id3 --> Jinjaa