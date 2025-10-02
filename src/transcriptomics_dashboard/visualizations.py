"""Visualization functions for transcriptomics data."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from typing import Optional, List, Tuple
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist


def create_volcano_plot(
    results: pd.DataFrame,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
    top_n_labels: int = 10,
    highlight_genes: Optional[List[str]] = None,
    title: str = "Volcano Plot"
) -> go.Figure:
    """
    Create interactive volcano plot.
    
    Args:
        results: DESeq2 results DataFrame
        fdr_threshold: FDR cutoff for significance
        lfc_threshold: Log2 fold change threshold
        top_n_labels: Number of top genes to label
        highlight_genes: List of specific genes to highlight
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    # Filter out genes with missing p-values
    plot_data = results.dropna(subset=['padj', 'log2FoldChange'])
    
    # Calculate -log10(padj)
    plot_data = plot_data.copy()
    plot_data['-log10padj'] = -np.log10(plot_data['padj'])
    
    # Replace infinite values
    max_log10p = plot_data['-log10padj'].replace([np.inf, -np.inf], np.nan).max()
    plot_data['-log10padj'] = plot_data['-log10padj'].replace([np.inf], max_log10p * 1.1)
    
    # Assign colors based on significance and direction
    colors = []
    for _, row in plot_data.iterrows():
        if row['padj'] < fdr_threshold and abs(row['log2FoldChange']) > lfc_threshold:
            if row['log2FoldChange'] > 0:
                colors.append('up')
            else:
                colors.append('down')
        else:
            colors.append('not_sig')
    
    plot_data['color'] = colors
    
    # Create color mapping
    color_map = {
        'up': '#E74C3C',      # Red
        'down': '#3498DB',    # Blue
        'not_sig': '#95A5A6'  # Gray
    }
    
    # Create figure
    fig = go.Figure()
    
    # Plot each category separately for better control
    for category, color in color_map.items():
        mask = plot_data['color'] == category
        data_subset = plot_data[mask]
        
        # Size by baseMean (expression level)
        sizes = np.log10(data_subset['baseMean'] + 1)
        sizes = (sizes / sizes.max() * 10) + 3  # Scale to 3-13
        
        fig.add_trace(go.Scatter(
            x=data_subset['log2FoldChange'],
            y=data_subset['-log10padj'],
            mode='markers',
            name=category.replace('_', ' ').title(),
            marker=dict(
                color=color,
                size=sizes,
                opacity=0.6 if category == 'not_sig' else 0.8,
                line=dict(width=0)
            ),
            text=data_subset['gene'],
            customdata=data_subset[['baseMean', 'pvalue', 'padj']],
            hovertemplate=(
                '<b>%{text}</b><br>' +
                'log2FC: %{x:.2f}<br>' +
                '-log10(padj): %{y:.2f}<br>' +
                'Padj: %{customdata[2]:.2e}<br>' +
                'BaseMean: %{customdata[0]:.0f}<br>' +
                '<extra></extra>'
            )
        ))
    
    # Add threshold lines
    fig.add_hline(
        y=-np.log10(fdr_threshold),
        line_dash="dash",
        line_color="gray",
        annotation_text=f"FDR = {fdr_threshold}",
        annotation_position="right"
    )
    
    fig.add_vline(
        x=lfc_threshold,
        line_dash="dash",
        line_color="gray"
    )
    
    fig.add_vline(
        x=-lfc_threshold,
        line_dash="dash",
        line_color="gray"
    )
    
    # Add labels for top genes
    if top_n_labels > 0:
        sig_genes = plot_data[plot_data['color'] != 'not_sig'].copy()
        sig_genes = sig_genes.sort_values('-log10padj', ascending=False)
        
        # Get top up and down regulated
        up_genes = sig_genes[sig_genes['color'] == 'up'].head(top_n_labels)
        down_genes = sig_genes[sig_genes['color'] == 'down'].head(top_n_labels)
        top_genes = pd.concat([up_genes, down_genes])
        
        for _, gene in top_genes.iterrows():
            fig.add_annotation(
                x=gene['log2FoldChange'],
                y=gene['-log10padj'],
                text=gene['gene'],
                showarrow=True,
                arrowhead=2,
                arrowsize=1,
                arrowwidth=1,
                arrowcolor='black',
                ax=20 if gene['log2FoldChange'] > 0 else -20,
                ay=-20,
                font=dict(size=9),
                bgcolor='rgba(255, 255, 255, 0.8)',
                borderpad=2
            )
    
    # Update layout
    fig.update_layout(
        title=title,
        xaxis_title="log<sub>2</sub> Fold Change",
        yaxis_title="-log<sub>10</sub> (adjusted p-value)",
        hovermode='closest',
        template='plotly_white',
        width=900,
        height=600,
        showlegend=True,
        legend=dict(
            x=0.02,
            y=0.98,
            bgcolor='rgba(255, 255, 255, 0.8)',
            bordercolor='black',
            borderwidth=1
        )
    )
    
    return fig


def create_heatmap(
    vst_counts: pd.DataFrame,
    de_results: pd.DataFrame,
    metadata: pd.DataFrame,
    top_n: int = 50,
    fdr_threshold: float = 0.05,
    cluster_genes: bool = True,
    cluster_samples: bool = True,
    title: str = "Expression Heatmap"
) -> go.Figure:
    """
    Create expression heatmap of top differentially expressed genes.
    
    Args:
        vst_counts: VST-transformed counts
        de_results: DESeq2 results
        metadata: Sample metadata
        top_n: Number of top genes to display
        fdr_threshold: FDR threshold for filtering
        cluster_genes: Whether to cluster genes
        cluster_samples: Whether to cluster samples
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    # Get significant genes
    sig_genes = de_results[de_results['padj'] < fdr_threshold].copy()
    sig_genes = sig_genes.sort_values('padj')
    
    # Select top N genes
    top_genes = sig_genes.head(top_n)['gene'].tolist()
    
    if len(top_genes) == 0:
        # No significant genes, use top by p-value
        top_genes = de_results.sort_values('pvalue').head(top_n)['gene'].tolist()
    
    # Subset VST counts
    heatmap_data = vst_counts.loc[top_genes].copy()
    
    # Z-score normalize
    heatmap_data = (heatmap_data.T - heatmap_data.mean(axis=1)) / heatmap_data.std(axis=1)
    heatmap_data = heatmap_data.T
    
    # Clustering
    gene_order = list(range(len(top_genes)))
    sample_order = list(range(len(heatmap_data.columns)))
    
    if cluster_genes and len(top_genes) > 1:
        gene_linkage = linkage(heatmap_data.values, method='average', metric='correlation')
        gene_dendro = dendrogram(gene_linkage, no_plot=True)
        gene_order = gene_dendro['leaves']
    
    if cluster_samples and len(heatmap_data.columns) > 1:
        sample_linkage = linkage(heatmap_data.T.values, method='average', metric='correlation')
        sample_dendro = dendrogram(sample_linkage, no_plot=True)
        sample_order = sample_dendro['leaves']
    
    # Reorder data
    heatmap_data = heatmap_data.iloc[gene_order, sample_order]
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data.values,
        x=heatmap_data.columns,
        y=heatmap_data.index,
        colorscale='RdBu_r',
        zmid=0,
        colorbar=dict(title="Z-score"),
        hovertemplate='Gene: %{y}<br>Sample: %{x}<br>Z-score: %{z:.2f}<extra></extra>'
    ))
    
    # Update layout
    fig.update_layout(
        title=title,
        xaxis_title="Samples",
        yaxis_title="Genes",
        template='plotly_white',
        width=800,
        height=max(400, len(top_genes) * 10),
        xaxis=dict(tickangle=-45),
        yaxis=dict(tickfont=dict(size=8))
    )
    
    return fig


def create_ma_plot(
    results: pd.DataFrame,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
    title: str = "MA Plot"
) -> go.Figure:
    """
    Create MA plot (mean expression vs log2 fold change).
    
    Args:
        results: DESeq2 results DataFrame
        fdr_threshold: FDR cutoff for significance
        lfc_threshold: Log2 fold change threshold
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    # Filter out genes with missing values
    plot_data = results.dropna(subset=['padj', 'log2FoldChange', 'baseMean']).copy()
    
    # Calculate log10(baseMean)
    plot_data['log10baseMean'] = np.log10(plot_data['baseMean'] + 1)
    
    # Assign colors
    colors = []
    for _, row in plot_data.iterrows():
        if row['padj'] < fdr_threshold and abs(row['log2FoldChange']) > lfc_threshold:
            if row['log2FoldChange'] > 0:
                colors.append('up')
            else:
                colors.append('down')
        else:
            colors.append('not_sig')
    
    plot_data['color'] = colors
    
    color_map = {
        'up': '#E74C3C',
        'down': '#3498DB',
        'not_sig': '#95A5A6'
    }
    
    # Create figure
    fig = go.Figure()
    
    for category, color in color_map.items():
        mask = plot_data['color'] == category
        data_subset = plot_data[mask]
        
        fig.add_trace(go.Scatter(
            x=data_subset['log10baseMean'],
            y=data_subset['log2FoldChange'],
            mode='markers',
            name=category.replace('_', ' ').title(),
            marker=dict(
                color=color,
                size=4,
                opacity=0.5 if category == 'not_sig' else 0.7,
                line=dict(width=0)
            ),
            text=data_subset['gene'],
            customdata=data_subset[['baseMean', 'padj']],
            hovertemplate=(
                '<b>%{text}</b><br>' +
                'log10(baseMean): %{x:.2f}<br>' +
                'log2FC: %{y:.2f}<br>' +
                'Padj: %{customdata[1]:.2e}<br>' +
                '<extra></extra>'
            )
        ))
    
    # Add threshold lines
    fig.add_hline(y=lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_hline(y=-lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_hline(y=0, line_color="black", line_width=1)
    
    # Update layout
    fig.update_layout(
        title=title,
        xaxis_title="log<sub>10</sub> (Mean Expression)",
        yaxis_title="log<sub>2</sub> Fold Change",
        hovermode='closest',
        template='plotly_white',
        width=900,
        height=600,
        showlegend=True
    )
    
    return fig


def create_pca_plot(
    vst_counts: pd.DataFrame,
    metadata: pd.DataFrame,
    condition_col: str,
    title: str = "PCA Plot"
) -> go.Figure:
    """
    Create PCA plot of samples.
    
    Args:
        vst_counts: VST-transformed counts
        metadata: Sample metadata
        condition_col: Column for coloring samples
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    from sklearn.decomposition import PCA
    
    # Transpose (samples as rows)
    data = vst_counts.T
    
    # Remove genes with zero variance
    data = data.loc[:, data.var() > 0]
    
    # Run PCA
    pca = PCA(n_components=min(10, data.shape[0], data.shape[1]))
    pca_coords = pca.fit_transform(data)
    
    # Create DataFrame
    pca_df = pd.DataFrame(
        pca_coords[:, :2],
        index=data.index,
        columns=['PC1', 'PC2']
    )
    
    # Add metadata
    pca_df = pca_df.join(metadata[[condition_col]])
    
    # Explained variance
    var_exp = pca.explained_variance_ratio_ * 100
    
    # Create plot
    fig = px.scatter(
        pca_df,
        x='PC1',
        y='PC2',
        color=condition_col,
        text=pca_df.index,
        title=title,
        labels={
            'PC1': f'PC1 ({var_exp[0]:.1f}%)',
            'PC2': f'PC2 ({var_exp[1]:.1f}%)'
        }
    )
    
    fig.update_traces(
        marker=dict(size=12, line=dict(width=1, color='white')),
        textposition='top center'
    )
    
    fig.update_layout(
        template='plotly_white',
        width=800,
        height=600,
        showlegend=True
    )
    
    return fig


if __name__ == "__main__":
    # Example usage with sample data
    import numpy as np
    
    np.random.seed(42)
    
    # Create sample DE results
    n_genes = 1000
    results = pd.DataFrame({
        'gene': [f'Gene_{i}' for i in range(n_genes)],
        'baseMean': np.random.lognormal(5, 2, n_genes),
        'log2FoldChange': np.random.normal(0, 2, n_genes),
        'lfcSE': np.random.uniform(0.1, 0.5, n_genes),
        'pvalue': np.random.beta(0.5, 5, n_genes),
        'padj': np.random.beta(0.5, 5, n_genes)
    })
    
    # Create volcano plot
    fig = create_volcano_plot(results)
    fig.show()