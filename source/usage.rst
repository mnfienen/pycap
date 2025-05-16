=====
Usage
=====
# Overview
The `pycap-dss` software provides analytical solutions to evaluate both streamflow depletion and drawdown resulting from pumping in groundwater wells. Several solutions for both depletion and drawdown are included and the modular software design allows for additional solutions to be added relatively easily.  

In the sections below, an overview of the capabilities for 

# Individual Solutions


# Well object
    - allows for superposition in time making time series pumping
    - multiple responses
    - changing the methods and apportionment

# Analysis Project
    - an example wrapping it all up together

# maybe geoprocessing?



Start by importing `pycap`

.. code-block:: python

    import pycap


Or import a well functions

.. code-block:: python

    import pycap.wells as wo

    ddwn = wo.theis()

Or import a specific function.

.. code-block:: python

    from pycap.wells import glover

    depl = glover()

    
