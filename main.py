import numericalunits as nu
from tests.runner import run
from tests.test_resonance import run_resonance_test

search_domain_kx = [-0.05/nu.nm, 0.05/nu.nm, 0, 0.4/nu.nm]
path = '1d_bragg.yaml'

if __name__ == "__main__":
    #test_davis()
    run(path=path, search_domain_kx=search_domain_kx, num_to_visualize=2, visualisation_type=['H', 'S'])
    #run_resonance_test(path)