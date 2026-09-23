import pandas as pd

def map_annotations(vertices_df, meta_df, map_by='cell', cols_to_map='all'):
    """
    Map annotations from meta_df to vertices_df based on map_by column.
    
    Parameters
    ----------
    vertices_df : DataFrame
        Input DataFrame with vertices and an id column.
    meta_df : DataFrame
        Metadata DataFrame with annotations and an id column.
    map_by : str
        Column name to map by (e.g., 'cell').
        
    Returns
    -------
    DataFrame
        Merged DataFrame with annotations.
    """
    
    # Merge the two dataframes on the specified column
    if cols_to_map != 'all':
        cols_to_map = [map_by] + cols_to_map
        meta_df = meta_df[cols_to_map]
    
    # Perform the merge
    merged_df = pd.merge(vertices_df, meta_df, on=map_by, how='left')
    
    return merged_df