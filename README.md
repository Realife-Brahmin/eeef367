# eeef367
**Theory in brief:** \
Transformers can develop 'Partial Discharge' or 'PD' (gas/liquid molecules from the transformer's own constituent materials, such as oil, becoming highly charged) due to regular operation in high voltage and high current conditions. These PD(s) can interact electrically with the flowing current and electric/magnetic fields within the transformer, disturbing its regular operation. Over time they can cause permanent efficiency losses and damage to the transformer in which they develop. Thus any PD(s) should ideally be neutralized upon detection. But they can be hard to locate due to the 'closed black-box' construction of the transformer. Given the terminal voltage and current readings (called the 'Line' and 'Neutral' terminals) of a transformer, the attached set of codes uses a 'Transfer Function Based Algorithm' in order to detect the location of a single PD within it, using the assumption that the transformer internal network can be represented as a homogeneous RLC Ladder network.

Code 1: Generalized Algorithm for generating Transfer Function of a given Transformer RLC ladder network.

Code 2: Detection of the location (node) of developed unknown Partial Discharge (modeled as time varying current) via Code 1.

This project was created as part of my undergraduate course **EEE F367: Lab Oriented Project** at **BITS Pilani Hyderabad Campus**. The algorithm is primarily inspired from [this paper](https://t.ly/xR6d "Research gate Link").
