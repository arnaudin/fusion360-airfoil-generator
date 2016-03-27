# fusion360-airfoil-generator

Generate sketch profiles for NACA airfoils in Fusion 360.

Available both as a script and add-in.

## Features :rocket:

- Supports 4 and 5 series NACA airfoils
- User specified number of segments
- Half cosine or linear spacing option
- Finite thickness trailing edge option

### Screenshots :camera: :computer:

![NACA Airfoil Generator](resources/app_screenshot_v1.jpg?raw=true "NACA Airfoil Generator")

![NACA Airfoil Rendering](resources/airfoil_rendering_v1.png?raw=true "NACA Airfoil Rendering")

## Installation and Usage :floppy_disk:

### Installation

Two options:

1. [Download](https://apps.autodesk.com/FUSION/en/Detail/Index?id=appstore.exchange.autodesk.com%3afusion360airfoilgenerator_windows64%3aen) and install from the Fusion 360 app store for Win64 only.
2. Download source files here, then see [how to install sample Add-Ins or Scripts](https://rawgit.com/AutodeskFusion360/AutodeskFusion360.github.io/master/Installation.html) within Fusion 360. Add-in version is recommended but script is also provided.  Be sure to copy the resources into a "resources" directory in the script or addin directory.  This process works with Mac OSX Fusion 360.

### Usage

1. If using the addin, select the *Airfoil* option from the *SKETCH* panel. For the script, execute script from *ADD-INS* panel.
2. Specify the NACA number for a 4 or 5 series airfoil.
3. Specify the number of points per side. Total points will be 2*numPoints+1. 
4. Selecting half cosine spacing will result in finer discretization near the leading edge of the airfoil compared to the default constant spacing. 
5. Selecting finite thickness will apply a finite thickness at the trailing edge compared to the default zero thickness.
6. Clicking *OK* generates the airfoil profile.
7. Scale, extrude, etc. in the modeling environment. 

## Credits :raised_hands:

NACA generation portions of this code are from [naca.py by Dirk Gorissen](https://github.com/dgorissen/naca). 

Thanks to the [Fusion 360 development team](https://github.com/AutodeskFusion360) for code samples and the development platform/API. 

## License :neckbeard:

This script is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT). Please see the [LICENSE](https://raw.githubusercontent.com/arnaudin/fusion360-airfoil-generator/master/LICENSE) file for full details.

## Author :pencil:
This add-in was created by [Ryan Arnaudin](http://ryanarnaudin.com).
