from tests.test_davis import test_davis
import numericalunits as nu
from tests.runner import run
from tests.test_reflection import run_reflection_test

search_domain_kx = [-0.05/nu.nm, 0.05/nu.nm, 0, 0.4/nu.nm]
path = 'air_metall_air.yaml'

if __name__ == "__main__":
    #test_davis()
    #run(path=path, num_to_visualize=5, visualisation_type=['H', 'S'])
    run_reflection_test(path)