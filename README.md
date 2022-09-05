# Kuramoto Neurons

This repository holds code related to "Synchronization in phase-coupled Kuramoto oscillator networks with axonal delay and synaptic plasticity" <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.032906>. It simulates emergent behavior of coupled oscillator networks with learning rules and time delays that more closely mimick the interactions between neurons in the brain than in traditional machine learning "neuronal network" based approaches.

![pic alt](./python/mygif.gif "polar plot showing synchronization of neurons into two different emergent groups")

## Igor Pro

It was originally written circa 2013 for my undergraduate honors thesis (available here: ) and the publication above. It as "procedure files" for a program called _Igor Pro_ which enabled easy and powerful visualization at the time but is unfortunately both locked-down/proprietary and obscure. This limited the potential use of the code and exploration of the models by others. I also didn't know how to use `git`, `python`, etc. at the time so this is just copied from surviving copies I had in an old backup. It is here just for the historical record.

## Rust

In 2022, I translated the Igor Pro code to the Rust programming language. Doing so provides 1. greater accessiblity of the model since rust is entirely open source and I can provide pre-compiled binaries, 2. much faster execution with modern optimizations and c-like performance, 3. a project to teach myself rust.
