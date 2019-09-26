from gnomad_hail import *
import pickle

public_root = 'gs://gnomad-public/papers/2019-high-ad-variants'
temp_root = 'gs://gnomad-tmp/high_ad_hets'


def build_frequency_case(mt):
    """
    Bin frequencies. Output refers to bin floor

    :param MatrixTable mt: Input MatrixTable with .freq annotation
    :return: Binned frequencies
    :rtype: FloatExpression
    """
    return (hl.case()
            .when(mt.freq[0].AF >= 0.2, 0.2)
            .when(mt.freq[0].AF >= 0.10, 0.10)
            .when(mt.freq[0].AF >= 0.05, 0.05)
            .when(mt.freq[0].AF >= 0.01, 0.01)
            .default(0))


def build_length_case(variant_length_expr):
    """
    Bin variant lengths: negative for deletions, positive for insertions, 0 for SNVs.
    Integer refers to bin floor

    :param IntExpression variant_length_expr: Input expression to bin
    :return: Binned variant lengths
    :rtype: IntExpression
    """
    return (hl.case()
            .when(variant_length_expr >= 10, 10)
            .when(variant_length_expr >= 5, 5)
            .when(variant_length_expr >= 1, 1)
            .when(variant_length_expr <= -10, -10)
            .when(variant_length_expr <= -5, -5)
            .when(variant_length_expr <= -1, -1)
            .default(0))


def format_hail_linreg(linreg_expr):
    """
    Format hail linear regression struct from arrays into slope and intercept

    :param StructExpression linreg_expr:
    :return: Formatted StructExpression
    :rtype: StructExpression
    """
    return linreg_expr.annotate(
        beta=linreg_expr.beta[1],
        standard_error=linreg_expr.standard_error[1],
        t_stat=linreg_expr.t_stat[1],
        p_value=linreg_expr.p_value[1],
        intercept_beta=linreg_expr.beta[0],
        intercept_standard_error=linreg_expr.standard_error[0],
        intercept_t_stat=linreg_expr.t_stat[0],
        intercept_p_value=linreg_expr.p_value[0]
    )


def main(args):
    hl.init(log='/high_ad_hets.log')

    for data_type in ('genomes', 'exomes'):
        if args.find_high_ad_hets:
            ht = get_gnomad_public_data(data_type)
            ht = ht.filter((ht.freq[0].AF >= 0.01) & (ht.freq[0].AF <= 0.5) & (hl.len(ht.filters) == 0))
            ht = ht.annotate(variant_length=hl.len(ht.alleles[1]) - hl.len(ht.alleles[0]))
            ht = ht.annotate(
                n_called_homs=ht.freq[0].homozygote_count,
                n_additional_homs=hl.sum(ht.ab_hist_alt.bin_freq[-2:]),  # corresponds to >= 0.9
                freq=build_frequency_case(ht),
                length_group=build_length_case(ht.variant_length)
            ).naive_coalesce(400)
            ht = ht.checkpoint(f'{public_root}/all_variants_{data_type}.ht', overwrite=args.overwrite, _read_if_exists=not args.overwrite)

            res = ht.group_by('freq', 'length_group').aggregate(
                called_homs=hl.agg.mean(ht.n_called_homs),
                added_homs=hl.agg.mean(ht.n_additional_homs),
                n_new_hom_variants=hl.agg.count_where((ht.n_called_homs == 0) & (ht.n_additional_homs > 0)),
                total_with_hom=hl.agg.count_where(ht.n_called_homs > 0),
                total_with_added_hom=hl.agg.count_where(ht.n_additional_homs > 0),
                total_with_1pct_added_hom=hl.agg.count_where(ht.n_additional_homs / ht.n_called_homs > 0.01),
                total_with_2pct_added_hom=hl.agg.count_where(ht.n_additional_homs / ht.n_called_homs > 0.02),
                total_with_5pct_added_hom=hl.agg.count_where(ht.n_additional_homs / ht.n_called_homs > 0.05),
                added_hom_hist=hl.agg.hist(ht.n_additional_homs, 0, 100, 100),
                prop_added_hom_hist=hl.agg.hist(hl.cond(ht.n_called_homs > 0, ht.n_additional_homs / ht.n_called_homs,
                                                        hl.int(ht.n_additional_homs > 0)), 0, 1, 100),
                total=hl.agg.count(),
                overall_ab_hist=hl.agg.array_agg(lambda ab: hl.agg.mean(ab), ht.ab_hist_alt.bin_freq)
            ).checkpoint(f'{public_root}/groupings_{data_type}.ht', overwrite=args.overwrite, _read_if_exists=not args.overwrite)
            res.flatten().export(f'{public_root}/groupings_{data_type}.txt')


        # Note: requires sample-level access to run
        if args.find_high_ad_hets_by_sample:
            mt = get_gnomad_data(data_type, non_refs_only=True, split=True, release_samples=True, release_annotations=True, adj=True, full_meta=True)
            mt = mt.filter_rows((hl.len(mt.filters) == 0) & (mt.freq[0].AF <= 0.5))
            mt = mt.annotate_rows(variant_length=hl.len(mt.alleles[1]) - hl.len(mt.alleles[0]),
                                  freq=build_frequency_case(mt))
            mt = mt.annotate_rows(length_group=build_length_case(mt.variant_length))
            mt = mt.group_rows_by('freq', 'length_group', 'segdup').partition_hint(400).aggregate_rows(
                n_variants=hl.agg.count()
            ).aggregate(
                number_reclassified=hl.agg.count_where(mt.GT.is_het() & (mt.AD[1] / mt.DP > 0.9)),
                number_homvars_with_ref_support=hl.agg.count_where(mt.GT.is_hom_var() & (mt.AD[0] > 0))
            )
            mt = mt.checkpoint(f'{temp_root}/num_reclassified_by_sample_{data_type}.mt', overwrite=args.overwrite, _read_if_exists=not args.overwrite)
            ht = mt.annotate_rows(mean_number_homvars_with_ref_support=hl.agg.mean(mt.number_homvars_with_ref_support),
                                  mean_number_reclassified=hl.agg.mean(mt.number_reclassified),
                                  linreg_number_reclassified=format_hail_linreg(hl.agg.linreg(mt.meta.freemix, [1.0, mt.number_reclassified])),
                                  linreg_number_homvars_with_ref_support=format_hail_linreg(hl.agg.linreg(mt.meta.freemix, [1.0, mt.number_homvars_with_ref_support]))
                                  ).rows()
            ht.flatten().export(f'{public_root}/contamination_summary_{data_type}.txt')

            ht = mt.select_cols(
                number_reclassified=hl.agg.sum(mt.number_reclassified),
                number_homvars_with_ref_support=hl.agg.sum(mt.number_homvars_with_ref_support),
                freemix=mt.meta.freemix
            ).cols().key_by().drop('s')
            ht.flatten().export(f'{temp_root}/num_reclassified_by_sample_all_variants_{data_type}.txt.bgz')
            data = ht.aggregate(linreg_number_reclassified=format_hail_linreg(hl.agg.linreg(ht.freemix, [1.0, ht.number_reclassified])),
                                linreg_number_homvars_with_ref_support=format_hail_linreg(hl.agg.linreg(ht.freemix, [1.0, ht.number_homvars_with_ref_support])))
            pprint(dict(data))
            with hl.hadoop_open(f'{public_root}/contamination_linregs_all_variants_{data_type}.pckl', 'wb') as f:
                pickle.dump(data, f)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--find_high_ad_hets', help='Find high AD hets from public data', action='store_true')
    parser.add_argument('--find_high_ad_hets_by_sample', help='Find high AD hets from individual level data', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite all data', action='store_true')
    args = parser.parse_args()

    main(args)
