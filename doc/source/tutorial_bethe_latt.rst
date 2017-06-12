Single-band Bethe lattice
-------------------------

In this example, we study a special one-band Hubbard model, 
which has semi-circular noninteracting density of states (dos),
It corresponds to Bethe lattice with infinite coordination number.

There is predefined class,
which helps generating the energy mesh with uniform weight.

.. autoclass:: pyglib.model.special.semicircular
    :members: 

In the model, we use half-band width as the energy unit. 
The noninteracting dos and cumulative dos is shown as below:

.. image:: _images/semicir_dos.png
    :alt: semicircular dos and cdos
    :scale: 100 %
    :align: center

A function to setup the model for *CyGutz* calculation 
has been defined,

.. autofunction:: pyglib.model.semicir.gutz_model_setup
    

