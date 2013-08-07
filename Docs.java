/**
 * @defgroup jrb JRB
 * @ingroup jarmosbase
 * @short Package for simulation of rbmatlab models.
 * 
 * This package is an extension/cleanup/rewrite of all the classes previously contained in the @ref rbappmit android
 * application. Even though basically everything has been structured and streamlined, the original work has been an
 * invaluable starting point in order to get grips on both Java RB simulations and android application development.
 * Hence, we owe a big thanks to David J. Knezevic and Phuong Huynh (MIT at that time) for the initial programming and
 * their agreement to use and improve the existing program.
 * 
 * Most of the included models stem from the old @ref rbappmit package, which is why there is an old model type "rbappmit"
 * for compatibility and a new type "JRB" which is to use for newly included models.
 * Check out the @ref jarmos_models section for details on how to include new models.
 * 
 * Check out the rbmatlab release on http://www.morepas.org/software/ for more information on how to generate and export
 * rbmatlab models to JaRMoS.
 * 
 * The JRB framework as a whole is published under the GNU GPL license stated below.
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>
 * 
 * @section rbappmit rbAppMIT
 * Some parts of this software have been taken from the "rbappmit" android application.
 * 
 * The author does not claim any copyright or originality of code and instead refers to the license conditions rbappmit
 * was published:
 * 
 * rbAPPmit: An Android front-end for the Certified Reduced Basis Method Copyright (C) 2010 David J. Knezevic and Phuong
 * Huynh
 * 
 * @ref rbappmit is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * @ref rbappmit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with rbAPPmit. If not, see
 * <http://www.gnu.org/licenses/>.
 * 
 * @{
 * 
 * @package rb
 * @short JRB main package
 * 
 * @package rb.affinefcn
 * @short Contains interfaces for affine functions
 * 
 * @package rb.test
 * @short Test classes
 * 
 * @package models.rbm_advec
 * @short Classes for time-dependent advection/diffusion problem
 * 
 * @package models.rbm_advec_cheatbase
 * @short Classes for time-dependent advection/diffusion problem with specially fitted basis
 * 
 * @package models.rbm_advec_tc
 * @short Classes for time-independent advection/diffusion problem
 * 
 * @}
 * 
 */
