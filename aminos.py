from scipy.stats import norm
import numpy as np 
import matplotlib.pyplot as plt

class Amino:
    
    def __init__(self, amino_label):

        self.amino_label = amino_label
        self.atoms_list = []
        self.gauss_mean = None
        self.gauss_std = None
    

    def make_gaussian_fit(self, slice_proportion = 0.9):

        cutoff_low = int(((1 - slice_proportion) / 2) * len(self.atoms_list))
        cutoff_high = len(self.atoms_list) - cutoff_low
        shift_val_list = [atom.shift_val for atom in self.atoms_list]
        self.gauss_mean, self.gauss_std = norm.fit(
            shift_val_list[cutoff_low:cutoff_high]
        )

    def find_outliers(self, num_sig_low, num_sig_high):

        for atom in self.atoms_list:
            diff_shift = abs(atom.shift_val - self.gauss_mean) 
            if (
                diff_shift >= num_sig_low * self.gauss_std and 
                diff_shift <= num_sig_high * self.gauss_std
            ):
                atom.is_outlier = True
    
    def make_plots(self, num_sig_low, num_sig_high):

        shift_val_list = [atom.shift_val for atom in self.atoms_list]
        fig, ax = plt.subplots()
        ax.hist(shift_val_list, 200, density=True)
        ax.axvspan(
            self.gauss_mean - num_sig_high * self.gauss_std, 
            self.gauss_mean - num_sig_low * self.gauss_std,
            color='yellow', alpha=0.5    
        )
        ax.axvspan(
            self.gauss_mean + num_sig_low * self.gauss_std, 
            self.gauss_mean + num_sig_high * self.gauss_std,
            color='yellow', alpha=0.5    
        )
        x = np.linspace(4, 12, 100)
        y = norm.pdf(x, self.gauss_mean, self.gauss_std)
        ax.plot(x, y)
        plt.show()
            

