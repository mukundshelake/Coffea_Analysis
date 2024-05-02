flowchart TD
    id1[(Data)] --> Skimmer --> id2{{skimmerOutput}} -->Extractor

    Extractor --> id3{{extractorOutput}}  -->Plotter
    %% Link *within* subgraph2:
    %% subgraph2 inherits the direction of the top-level graph (LR)
    unfilledTex --> Jinja --> filledReport

    id3 --> Jinjaa
