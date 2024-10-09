# Instructions for Configuring and Running Quicksilver

First, change the path of your Power API install in the [Makefile](./src/Makefile) to your install directory.

Then build with `make` and run `source xml_profile`.

Next, to run with RAPL, use the following command: `sudo --preserve-env ./src/qs`.
