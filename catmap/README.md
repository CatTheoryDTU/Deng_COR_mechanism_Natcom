Microkinetic model for estimating the H and CO coverage
Catmap version 0.3.1 was used cloned from ```git@github.com:sringe/catmap-1.git``` was used.

The base directory is ```catmap```. In order to run the microkinetic models simply run

```python3 analyze_mkm.py 100```

```python3 analyze_mkm.py 211```

The results will be written into ```../results```

Note that by changing :dattype: in analyze_mkm.py the TOF's can be plotted as
well. Furthermore, by adding ```'CO_'+facet``` in the cov_and_rate dict, the CO
coverage will be plotted as well.
Careful! If the pkl file ```data.pkl``` is not present in the respective facet
directories you might need to ramp the adsorbate-adsorbate interactions or
turn them off in :reactions.mkm:
