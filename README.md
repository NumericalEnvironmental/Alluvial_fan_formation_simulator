# Alluvial_fan_formation_simulator

![Preview](https://numericalenvironmental.files.wordpress.com/2017/12/gravel_surface.png?w=1158&h=494)

This repository contains the source code, example input files, and Windows executable for an alluvial fan formation simulator, written in D language. Some of the technical elements and assumptions used in the model are discussed in a post on my blog, https://numericalenvironmental.wordpress.com/2017/12/28/simulating-the-formation-of-an-alluvial-fan/. Here is a brief summary of the code itself.

The file alluvium.d contains several classes that handle different facets of the overall algorithm. These include:

* Sediment class - use to specify the properties of individual sediment size class components (e.g., clay-size particles, sand-size particles, etc.) which include size distributions (d50 and d90) and a method for estimating settling velocity
* Cell class - spatial element within the model domain, characterized by a location, bottom (initial) elevation, and a dynamic stack of sediments consisting of slices of uniform thickness, with each containing its own blend of sediment size classes
* Stream class - methods to address flow hydraulics, sediment load, deposition rates along a reach
* Grid class - methods to implement a meandering random walk algorithm, and for assigning cells to a stream/reach
* Drainage class - sediment distribution (location, influx, sediment composition and variability)
* Discharge class - runoff generator; intensity and duration

The separate utility.d file contains a few short miscellaneous standalone functions (not associated with any class) used to extract various kinds of information from arrays.

Please note that the Windows executable was compiled using dmd version 2.077. The output of the model was distinctly different and incorrect when compiled under an earlier version (for reasons I could not determine and did not pursue), so be sure to use an up-to-date compiler release if you choose to modify and recompile on your own.

Required input (text) files to run the code include the following:

* Cells.txt - cell numbering scheme and associated starting (pre-sedimentation) elevations
* Discharge.txt - lognormal distribution constraints for discharge event intensity and duration
* Drainage.txt - drainage locations, plus constraints used to model sediment composition captured in influx onto fan
* Grid_params.txt - grid size and other global parameters (e.g., weighting factor for discharge momentum) 
* Sediments.txt - sediment size class parameters (e.g., d50, d90)
* Stream.txt - stream attributes (Manning roughness coefficient, constraints on changes in cell sediment thickness during a given time step, etc.)

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

