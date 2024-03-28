import pandas as pd

def calculate_combined_variance(df, standard_deviation = False):
    """
    Calculate the combined standard deviation from multiple groups.
    
    Parameters:
    - df: pandas DataFrame with columns 'mean', 'variance', and 'sample_size'.
    
    Returns:
    - combined_std: Combined standard deviation of the groups.
    """
    # Calculate the weighted average of the means (grand mean)
    grand_mean = (df['mean'] * df['sample_size']).sum() / df['sample_size'].sum()
    weighted_sum_of_variances = ((df['variance'] * (df['sample_size'] - 1)).sum() +
                                 (df['sample_size'] * ((df['mean'] - grand_mean) ** 2)).sum())

    # Calculate the combined variance or standard deviation
    if standard_deviation:
        print('returning combined standard deviation')
        combined_std_var = combined_variance ** 0.5
    
    else:
        print('returning combined variance')
        combined_std_var = weighted_sum_of_variances / (df['sample_size'].sum() - len(df))

    return combined_std
