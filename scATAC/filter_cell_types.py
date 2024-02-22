import click
import pandas as pd

@click.command
@click.option("--cell_frag_counts", required=True)
@click.option("--min_frags", default=2e6)
@click.option("--output", required=True)
def main(cell_frag_counts, min_frags, output):
	df = pd.read_csv(cell_frag_counts, sep='\t')
	filtered_df = df[df['num_fragments'] >= min_frags]
	filtered_df.to_csv(output, sep='\t', index=False)

if __name__ == "__main__":
	main()