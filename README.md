
<h1 align="center">fvhyper</h1>
<p align="center"><i>Finite volume framework for high performance</i></p>
<br>

<!-- Badges -->
<p align="center">
  <a href="https://github.com/Alexsimulation/fvhyper/graphs/contributors">
    <img src="https://img.shields.io/github/contributors/Alexsimulation/fvhyper" alt="contributors" />
  </a>
  <a href="">
    <img src="https://img.shields.io/github/last-commit/Alexsimulation/fvhyper" alt="last update" />
  </a>
  <a href="https://github.com/Alexsimulation/fvhyper/network/members">
    <img src="https://img.shields.io/github/forks/Alexsimulation/fvhyper" alt="forks" />
  </a>
  <a href="https://github.com/Alexsimulation/fvhyper/stargazers">
    <img src="https://img.shields.io/github/stars/Alexsimulation/fvhyper" alt="stars" />
  </a>
  <a href="https://github.com/Alexsimulation/fvhyper/issues/">
    <img src="https://img.shields.io/github/issues/Alexsimulation/fvhyper" alt="open issues" />
  </a>
  <a href="https://github.com/Alexsimulation/fvhyper/blob/master/LICENSE">
    <img src="https://img.shields.io/github/license/Alexsimulation/fvhyper.svg" alt="license" />
  </a>
</p>

<div align="center">
<h4>
    <a href="https://github.com/Alexsimulation/fvhyper/wiki">Documentation</a>
  <span> · </span>
    <a href="https://github.com/Alexsimulation/fvhyper/issues/">Report Bug</a>
  <span> · </span>
    <a href="https://github.com/Alexsimulation/fvhyper/issues/">Request Feature</a>
</h4>
</div>

<br>



## :star2: About the project

This project aims to provide a MPI parallelizable framework to build finite volume solvers. For now, the scope of this project is to provie 2D solver capabilities, but it might be extended to 3D someday.

In its current form, the user can specify one partial differential equation with any number of variables, and any number of boundary conditions. To enter the pde, the user must define the flux function in the finite volume formulation.

## :wrench: Installation

To install, clone the git repository into your installation directory of choice, and then export a variable FVHYPER_DIR with the path to this installation to your .bashrc. Per example, if I installed the files to a $HOME/softwares directory:

```
echo "export FVHYPER_DIR=$HOME/softwares/fvhyper" >> ~/.bashrc
```

Then, using *make*, a minimal makefile contents for a given script *main.cpp* should look like this:

```
SOURCES := $(shell find $(FVHYPER_SOURCEDIR)/fvhyper/src -name '*.cpp')

INCLUDES := -I${FVHYPER_DIR}

build:
	mpic++ -o main main.cpp ${SOURCES} ${INCLUDES}

```

INCLUDES := -I${FVHYPER_DIR}

## :light_rail: Roadmap

 * [x] Mesh reader
 * [x] Explicit Euler solver
 * [x] PVTU file export
 * [x] Example 2D tests
 * [ ] Explicit Rk4 solver
 * [ ] Implicit Euler solver
 * [ ] 3D solvers

## :hand: Contribute

Contributions are always welcome! Please create a PR to add Github Profile.

## :pencil: License

This project is licensed under [MIT](https://opensource.org/licenses/MIT) license.

## :man_astronaut: Show your support

Give a ⭐️ if this project helped you!
