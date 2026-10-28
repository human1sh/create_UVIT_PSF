import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt

# Patch Matplotlib cm.get_cmap for legacy libraries like aplpy
if not hasattr(mpl.cm, 'get_cmap'):
    mpl.cm.get_cmap = plt.get_cmap

def setup_environment():
    """Configures environment parameters and suppresses warnings."""
    # Suppress Gaia release notices and Astropy deprecation warnings via string matching
    warnings.filterwarnings("ignore", message=".*current archive release.*")
    warnings.filterwarnings("ignore", message=".*TestRunner.*")
    warnings.filterwarnings("ignore", message=".*astropy.tests.runner.*")
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    
    # Matplotlib styling configurations
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.size'] = 14
    mpl.rcParams["figure.figsize"] = (8, 6)
    mpl.rcParams['savefig.dpi'] = 300
    mpl.rcParams["legend.fontsize"] = 12
    mpl.rcParams["axes.grid"] = False
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True