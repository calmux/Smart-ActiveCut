! [logo] (activeCut_logo.png)\n## Smart-ActiveCut\nActive learning-based interactive tool for semi-supervised image segmentation (lesion\/tumor segmentation as an example of application)\n\n\n### Overview of the proposed algorithm\n! [flowchart] (activeCut_flowchart.png)\n\n\n### Install & usage:\n\nTo compile the code:\n\n1) install cmake, boost and ITK libraries on your machine\n\n2) create a new directory for out-of-source build\n```sh\n$ mkdir build\n$ cd build\n$ ccmake ..\/src\n```\n3) 'c' to configue and 'g' to generate makefile.\n\n4) make\n\nAfter step 4, just type \n```sh\n$ .\/activeCutSeg\n```\nInput arguments will be printed.\n\nExample command to run activeCut:\n```sh\n$ .\/activeCutSeg -d ..\/allcha