from tests.test_davis import test_davis
import numericalunits as nu
from tests.runner import run

search_domain_kx = [-0.05/nu.nm, 0.05/nu.nm, 0, 0.4/nu.nm]
path = 'davis.yaml'

if __name__ == "__main__":
    #test_davis()
    run(path=path, search_domain_kx=search_domain_kx, num_to_visualize=5)