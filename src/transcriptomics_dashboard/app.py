"""Main Dash application for transcriptomics dashboard."""

import base64
import io
import logging
from typing import Optional, Dict

import dash
from dash import dcc, html, Input, Output, State, callback_context
import dash_bootstrap_components as dbc
import pandas as pd

from transcriptomics_dashboard.config import get_config
from transcriptomics_dashboard.validation import (
    read_count_matrix,
    read_metadata,
    validate_analysis_inputs
)
from transcriptomics_dashboard.deseq2 import run_deseq2, DESeq2Error
from transcriptomics_dashboard.visualizations import (
    create_volcano_plot,
    create_heatmap,
    create_ma_plot,
    create_pca_plot
)


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize config
config = get_config()

# Initialize Dash app
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    title=config.app_title
)

# Store for session data
session_data = {}


def create_layout():
    """Create the main dashboard layout."""
    return dbc.Container([
        # Header
        dbc.Row([
            dbc.Col([
                html.H1(config.app_title, className="text-primary mb-2"),
                html.H4("Interactive RNA-seq Differential Expression Analysis", 
                       className="text-secondary mb-4"),
                html.Hr()
            ])
        ]),
        
        # Upload Section
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(html.H5("1. Upload Data", className="mb-0")),
                    dbc.CardBody([
                        html.Label("Count Matrix:", className="fw-bold"),
                        dcc.Upload(
                            id='upload-counts',
                            children=html.Div([
                                'Drag and Drop or ',
                                html.A('Select Count Matrix File')
                            ]),
                            style={
                                'width': '100%',
                                'height': '60px',
                                'lineHeight': '60px',
                                'borderWidth': '2px',
                                'borderStyle': 'dashed',
                                'borderRadius': '5px',
                                'textAlign': 'center',
                                'margin': '10px 0'
                            },
                            multiple=False
                        ),
                        html.Div(id='counts-upload-status', className="mt-2"),
                        
                        html.Hr(className="my-3"),
                        
                        html.Label("Metadata:", className="fw-bold mt-3"),
                        dcc.Upload(
                            id='upload-metadata',
                            children=html.Div([
                                'Drag and Drop or ',
                                html.A('Select Metadata File')
                            ]),
                            style={
                                'width': '100%',
                                'height': '60px',
                                'lineHeight': '60px',
                                'borderWidth': '2px',
                                'borderStyle': 'dashed',
                                'borderRadius': '5px',
                                'textAlign': 'center',
                                'margin': '10px 0'
                            },
                            multiple=False
                        ),
                        html.Div(id='metadata-upload-status', className="mt-2"),
                        
                        html.Hr(className="my-3"),
                        
                        dbc.Button(
                            "Download Example Data",
                            id="download-example-btn",
                            color="secondary",
                            size="sm",
                            className="mt-2"
                        )
                    ])
                ], className="mb-4")
            ], width=12)
        ]),
        
        # Parameters Section
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(html.H5("2. Analysis Parameters", className="mb-0")),
                    dbc.CardBody([
                        html.Label("Condition Column:", className="fw-bold"),
                        dcc.Dropdown(
                            id='condition-column-dropdown',
                            placeholder="Select condition column from metadata",
                            disabled=True
                        ),
                        
                        html.Label("Comparison:", className="fw-bold mt-3"),
                        dbc.Row([
                            dbc.Col([
                                dcc.Dropdown(
                                    id='numerator-dropdown',
                                    placeholder="Numerator",
                                    disabled=True
                                )
                            ], width=5),
                            dbc.Col(html.Div("vs", className="text-center pt-2"), width=2),
                            dbc.Col([
                                dcc.Dropdown(
                                    id='denominator-dropdown',
                                    placeholder="Denominator",
                                    disabled=True
                                )
                            ], width=5)
                        ]),
                        
                        html.Hr(className="my-3"),
                        
                        dbc.Row([
                            dbc.Col([
                                html.Label("FDR Threshold:", className="fw-bold"),
                                dbc.Input(
                                    id='fdr-input',
                                    type='number',
                                    value=config.defaults.fdr_threshold,
                                    min=0,
                                    max=1,
                                    step=0.01
                                )
                            ], width=6),
                            dbc.Col([
                                html.Label("Log2FC Threshold:", className="fw-bold"),
                                dbc.Input(
                                    id='lfc-input',
                                    type='number',
                                    value=config.defaults.log2fc_threshold,
                                    min=0,
                                    step=0.1
                                )
                            ], width=6)
                        ]),
                        
                        html.Hr(className="my-3"),
                        
                        dbc.Button(
                            "Run Analysis",
                            id="run-analysis-btn",
                            color="primary",
                            size="lg",
                            className="w-100 mt-3",
                            disabled=True
                        ),
                        
                        html.Div(id='analysis-status', className="mt-3")
                    ])
                ], className="mb-4")
            ], width=12)
        ]),
        
        # Results Section
        dbc.Row([
            dbc.Col([
                dcc.Loading(
                    id="loading-results",
                    type="default",
                    children=html.Div(id='results-container')
                )
            ])
        ]),
        
        # Hidden divs for storing data
        html.Div(id='counts-data-store', style={'display': 'none'}),
        html.Div(id='metadata-data-store', style={'display': 'none'}),
        html.Div(id='analysis-results-store', style={'display': 'none'})
        
    ], fluid=True, className="py-4")


app.layout = create_layout()


# Callbacks

@app.callback(
    [Output('counts-upload-status', 'children'),
     Output('counts-data-store', 'children')],
    Input('upload-counts', 'contents'),
    State('upload-counts', 'filename')
)
def upload_counts(contents, filename):
    """Handle count matrix upload."""
    if contents is None:
        return "", ""
    
    try:
        # Parse uploaded file
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        
        # Read file
        if filename.endswith('.csv'):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), index_col=0)
        elif filename.endswith('.tsv') or filename.endswith('.txt'):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t', index_col=0)
        elif filename.endswith(('.xlsx', '.xls')):
            df = pd.read_excel(io.BytesIO(decoded), index_col=0)
        else:
            return dbc.Alert("Unsupported file format", color="danger"), ""
        
        # Validate
        from transcriptomics_dashboard.validation import validate_count_matrix
        result, schema = validate_count_matrix(df)
        
        # Store data
        session_data['counts'] = df
        
        # Create status message
        if result.valid:
            msg = dbc.Alert([
                html.Strong("✓ Count matrix uploaded successfully!"),
                html.Br(),
                f"Genes: {schema.n_genes:,} | Samples: {schema.n_samples}"
            ], color="success")
        else:
            msg = dbc.Alert([
                html.Strong("✗ Validation errors:"),
                html.Ul([html.Li(err) for err in result.errors])
            ], color="danger")
        
        if result.warnings:
            msg = html.Div([
                msg,
                dbc.Alert([
                    html.Strong("⚠ Warnings:"),
                    html.Ul([html.Li(w.message) for w in result.warnings])
                ], color="warning")
            ])
        
        return msg, "loaded"
        
    except Exception as e:
        logger.error(f"Error uploading counts: {e}")
        return dbc.Alert(f"Error: {str(e)}", color="danger"), ""


@app.callback(
    [Output('metadata-upload-status', 'children'),
     Output('metadata-data-store', 'children'),
     Output('condition-column-dropdown', 'options'),
     Output('condition-column-dropdown', 'disabled')],
    Input('upload-metadata', 'contents'),
    State('upload-metadata', 'filename')
)
def upload_metadata(contents, filename):
    """Handle metadata upload."""
    if contents is None:
        return "", "", [], True
    
    try:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        
        if filename.endswith('.csv'):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), index_col=0)
        elif filename.endswith('.tsv') or filename.endswith('.txt'):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t', index_col=0)
        elif filename.endswith(('.xlsx', '.xls')):
            df = pd.read_excel(io.BytesIO(decoded), index_col=0)
        else:
            return dbc.Alert("Unsupported file format", color="danger"), "", [], True
        
        # Store data
        session_data['metadata'] = df
        
        # Get column options
        options = [{'label': col, 'value': col} for col in df.columns]
        
        msg = dbc.Alert([
            html.Strong("✓ Metadata uploaded successfully!"),
            html.Br(),
            f"Samples: {len(df)} | Columns: {len(df.columns)}"
        ], color="success")
        
        return msg, "loaded", options, False
        
    except Exception as e:
        logger.error(f"Error uploading metadata: {e}")
        return dbc.Alert(f"Error: {str(e)}", color="danger"), "", [], True


@app.callback(
    [Output('numerator-dropdown', 'options'),
     Output('denominator-dropdown', 'options'),
     Output('numerator-dropdown', 'disabled'),
     Output('denominator-dropdown', 'disabled'),
     Output('run-analysis-btn', 'disabled')],
    Input('condition-column-dropdown', 'value'),
    State('metadata-data-store', 'children'),
    State('counts-data-store', 'children')
)
def update_condition_options(condition_col, metadata_loaded, counts_loaded):
    """Update condition value dropdowns."""
    if not condition_col or not metadata_loaded or not counts_loaded:
        return [], [], True, True, True
    
    metadata = session_data.get('metadata')
    if metadata is None or condition_col not in metadata.columns:
        return [], [], True, True, True
    
    # Get unique conditions
    conditions = sorted(metadata[condition_col].unique())
    options = [{'label': c, 'value': c} for c in conditions]
    
    return options, options, False, False, False


@app.callback(
    [Output('analysis-status', 'children'),
     Output('analysis-results-store', 'children'),
     Output('results-container', 'children')],
    Input('run-analysis-btn', 'n_clicks'),
    State('condition-column-dropdown', 'value'),
    State('numerator-dropdown', 'value'),
    State('denominator-dropdown', 'value'),
    State('fdr-input', 'value'),
    State('lfc-input', 'value')
)
def run_analysis(n_clicks, condition_col, numerator, denominator, fdr, lfc):
    """Run DESeq2 analysis."""
    if not n_clicks:
        return "", "", ""
    
    if not all([condition_col, numerator, denominator]):
        return dbc.Alert("Please select all analysis parameters", color="warning"), "", ""
    
    try:
        counts = session_data.get('counts')
        metadata = session_data.get('metadata')
        
        if counts is None or metadata is None:
            return dbc.Alert("Please upload both count matrix and metadata", color="danger"), "", ""
        
        # Validate inputs
        from transcriptomics_dashboard.validation import validate_analysis_inputs
        validation = validate_analysis_inputs(counts, metadata, condition_col)
        
        if not validation.valid:
            return dbc.Alert([
                html.Strong("Validation failed:"),
                html.Ul([html.Li(err) for err in validation.errors])
            ], color="danger"), "", ""
        
        # Run DESeq2
        logger.info("Starting DESeq2 analysis...")
        results = run_deseq2(
            counts,
            metadata,
            condition_col,
            contrast=(numerator, denominator),
            fdr_threshold=fdr,
            lfc_threshold=lfc
        )
        
        # Store results
        session_data['de_results'] = results['results']
        session_data['normalized_counts'] = results['normalized_counts']
        session_data['vst_counts'] = results['vst_counts']
        
        # Create visualizations
        volcano_fig = create_volcano_plot(
            results['results'],
            fdr_threshold=fdr,
            lfc_threshold=lfc,
            title=f"Volcano Plot: {numerator} vs {denominator}"
        )
        
        heatmap_fig = create_heatmap(
            results['vst_counts'],
            results['results'],
            metadata,
            top_n=50,
            fdr_threshold=fdr
        )
        
        ma_fig = create_ma_plot(
            results['results'],
            fdr_threshold=fdr,
            lfc_threshold=lfc
        )
        
        pca_fig = create_pca_plot(
            results['vst_counts'],
            metadata,
            condition_col
        )
        
        # Count significant genes
        sig_up = (results['results']['direction'] == 'up').sum()
        sig_down = (results['results']['direction'] == 'down').sum()
        
        # Create results layout
        results_layout = html.Div([
            dbc.Alert([
                html.H5("✓ Analysis Complete!", className="alert-heading"),
                html.Hr(),
                html.P([
                    f"Significant genes (FDR < {fdr}, |log2FC| > {lfc}): ",
                    html.Strong(f"{sig_up + sig_down}")
                ]),
                html.P([
                    f"Up-regulated: {sig_up} | Down-regulated: {sig_down}"
                ])
            ], color="success"),
            
            dbc.Tabs([
                dbc.Tab(label="Volcano Plot", children=[
                    dcc.Graph(figure=volcano_fig, style={'height': '600px'})
                ]),
                dbc.Tab(label="Heatmap", children=[
                    dcc.Graph(figure=heatmap_fig)
                ]),
                dbc.Tab(label="MA Plot", children=[
                    dcc.Graph(figure=ma_fig, style={'height': '600px'})
                ]),
                dbc.Tab(label="PCA", children=[
                    dcc.Graph(figure=pca_fig, style={'height': '600px'})
                ]),
                dbc.Tab(label="Results Table", children=[
                    dbc.Table.from_dataframe(
                        results['results'].head(100)[
                            ['gene', 'baseMean', 'log2FoldChange', 'pvalue', 'padj', 'direction']
                        ].round(4),
                        striped=True,
                        bordered=True,
                        hover=True,
                        size='sm',
                        className="mt-3"
                    ),
                    dbc.Button(
                        "Download Full Results (CSV)",
                        id="download-results-btn",
                        color="primary",
                        className="mt-3"
                    ),
                    dcc.Download(id="download-results")
                ])
            ])
        ])
        
        status = dbc.Alert("Analysis completed successfully!", color="success")
        
        return status, "complete", results_layout
        
    except DESeq2Error as e:
        logger.error(f"DESeq2 error: {e}")
        return dbc.Alert(f"DESeq2 Error: {str(e)}", color="danger"), "", ""
    except Exception as e:
        logger.error(f"Analysis error: {e}", exc_info=True)
        return dbc.Alert(f"Error: {str(e)}", color="danger"), "", ""


@app.callback(
    Output("download-results", "data"),
    Input("download-results-btn", "n_clicks"),
    prevent_initial_call=True
)
def download_results(n_clicks):
    """Download results as CSV."""
    if n_clicks and 'de_results' in session_data:
        return dcc.send_data_frame(
            session_data['de_results'].to_csv,
            "de_results.csv",
            index=False
        )


def main():
    """Run the Dash application."""
    app.run_server(
        debug=config.debug,
        host=config.host,
        port=config.port
    )


if __name__ == '__main__':
    main()