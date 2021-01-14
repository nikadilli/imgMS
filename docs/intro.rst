Introduction
============

``imgMS`` is a Python package which aims to provide an easy and intuitive way of working with LA-ICP-MS data acquired into single file. 

Motivation
**********

As a Phd student in the fiels of analytical chemistry using LA-ICP-MS for rutine analysis as well as research tasks, I found it pretty hard and time consuming to evaluate data without a appropriate software. Despite the fact that there exists multiple free or open source software, I couldn't find any that would fit our workflow. What started as an simple script for data reduction of LA-ICP-MS data is now a complete python package.  

This package is intended to provide a quick, as well as (hopefully) easy to undestand, way of getting a results from bulk analysis or elemental imaging using LA-ICP-MS.

Limitations
***********

- So far there is only option to directly import data from 2 instruments (ThermoFisher Element2 and Agilent) in 3 formats (.csv, .xlsx, .asc). In the future there will be more options.

- imgMS only works with multiple analysis with reference analysis acquired in one file, unlike many other software aimed on LA-ICP-MS data reduction. 

