import torch

# One-hot encoding
def onehot(sequence: str):
    """
    Code from Handling variable size DNA inputs
    https://elferachid.medium.com/handling-variable-size-dna-input-c4bbb1df0458
    """
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    encoded_seq = [mapping[nt] for nt in sequence]
    # print(np.eye(4)[encoded_seq].shape)
    return torch.eye(4)[encoded_seq]


def one_hot_encode_sequence(sequence):
    """
    Given a DNA sequence, return its one-hot encoding
    """
    # Dictionary returning one-hot encoding for each nucleotide
    nuc_d = {'dA':[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
             'dCm':[0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
             'dG':[0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
             'dT':[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
             'cEtA':[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
             'cEtCm':[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
             'cEtG':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
             'cEtT':[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]}
    # Create array from nucleotide sequence
    return torch.tensor([nuc_d[x] for x in sequence])



# import category_encoders