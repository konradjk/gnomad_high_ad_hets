from gnomad_hail import *
import pickle

public_root = 'gs://gnomad-public/papers/2019-high-ad-variants'
temp_root = 'gs://gnomad-tmp/high_ad_hets'
POPS_TO_USE = ['fin', 'nwe', 'onf', 'swe', 'est']


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

    if args.recompute_ccr5_homozygosity:
        hts = []
        for data_type in ('genomes', 'exomes'):
            mt = get_gnomad_data(data_type, adj=True, release_samples=True, raw=True, split=False, full_meta=True)
            mt = mt.filter_rows(mt.locus == hl.parse_locus('3:46414943', 'GRCh37')).select_entries('GT', 'AD', 'DP', 'GQ', 'PL')
            mt = mt.checkpoint(f'{temp_root}/ccr5_{data_type}.mt', overwrite=args.overwrite, _read_if_exists=not args.overwrite)
            mt = hl.filter_alleles_hts(mt, lambda allele, _: hl.is_indel(mt.alleles[0], allele), subset=True)
            mt = annotate_adj(mt)
            if args.remove_contaminated:
                mt = mt.filter_cols(mt.meta.freemix < 0.01)
            mt = mt.select_cols(pop=hl.or_else(mt.meta.subpop, mt.meta.pop)).select_rows()
            hts.append(mt.entries())

        ht = hts[0].union(hts[1])
        ht = ht.filter(hl.array(POPS_TO_USE).contains(ht.pop))
        ht = ht.union(ht.annotate(pop='all'))

        annotations = {
            'reported': hl.agg.call_stats(ht.GT, ['', '']),
            'actual': hl.agg.call_stats(hl.cond(ht.GT.is_het() & (ht.AD[1] / ht.DP > 0.9), hl.call(1, 1), ht.GT), ['', ''])
        }
        ht = ht.group_by(ht.pop).aggregate(**annotations)
        ht = ht.annotate(
            **{f'{ann}_expected_homs': ht[ann].AN * (ht[ann].AF[1] ** 2) / 2 for ann in annotations.keys()}
        )
        ht = ht.annotate(
            **{f'{ann}_bi': ht[ann].homozygote_count[1] / ht[f'{ann}_expected_homs'] for ann in annotations.keys()}
        )
        fname = f'{temp_root}/ccr5_results{"_nocontam" if args.remove_contaminated else ""}.ht'
        ht = ht.checkpoint(fname, overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        ht.select(reported_homs=ht.reported.homozygote_count[1], actual_homs=ht.actual.homozygote_count[1],
                  reported_AF=ht.reported.AF[1], actual_AF=ht.actual.AF[1],
                  reported_expected_homs=ht.reported_expected_homs, actual_expected_homs=ht.actual_expected_homs,
                  reported_bi=ht.reported_bi, actual_bi=ht.actual_bi).show(20)

    if args.recompute_background_distribution:
        POPS_TO_USE.append('all')
        pop_ht = hl.read_table(f'{temp_root}/ccr5_results.ht')
        ccr5_freqs = pop_ht.aggregate(hl.dict(hl.agg.collect((pop_ht.pop, pop_ht.row_value))), _localize=False)

        exome_mt = get_gnomad_data('exomes', adj=True, non_refs_only=True, release_samples=True, release_annotations=True)

        def format_pop_name(pop):
            if pop == 'all':
                return 'gnomad'
            elif pop == 'fin':
                return 'gnomad_fin'
            else:
                return 'gnomad_nfe_' + pop

        exome_mt = exome_mt.filter_rows(hl.any(lambda x: x >= 0.05,
            [exome_mt.freq[exome_mt.freq_index_dict[format_pop_name(pop)]].AF for pop in POPS_TO_USE]
        ))
        exome_mt = exome_mt.select_cols(pop=hl.or_else(exome_mt.meta.subpop, exome_mt.meta.pop))
        annotations = {
            'reported': hl.agg.call_stats(exome_mt.GT, 2),
            'actual': hl.agg.call_stats(hl.cond(exome_mt.GT.is_het() & (exome_mt.AD[1] / exome_mt.DP > 0.9),
                                                hl.call(1, 1), exome_mt.GT), 2)
        }
        exome_mt = exome_mt.group_cols_by(exome_mt.pop).aggregate(**annotations)
        exome_mt = exome_mt.checkpoint(f'{temp_root}/exome_common_freqs.mt', overwrite=args.overwrite, _read_if_exists=not args.overwrite)

        genome_mt = get_gnomad_data('genomes', adj=True, non_refs_only=True, release_samples=True, release_annotations=True)
        genome_mt = genome_mt.filter_rows(hl.is_defined(exome_mt.rows()[genome_mt.row_key]))
        genome_mt = genome_mt.select_cols(pop=hl.or_else(genome_mt.meta.subpop, genome_mt.meta.pop))
        annotations = {
            'reported': hl.agg.call_stats(genome_mt.GT, 2),
            'actual': hl.agg.call_stats(hl.cond(genome_mt.GT.is_het() & (genome_mt.AD[1] / genome_mt.DP > 0.9),
                                                hl.call(1, 1), genome_mt.GT), 2)
        }
        genome_mt = genome_mt.group_cols_by(genome_mt.pop).aggregate(**annotations)
        genome_mt = genome_mt.checkpoint(f'{temp_root}/genome_common_freqs.mt', overwrite=args.overwrite, _read_if_exists=not args.overwrite)

        genome_rows_join = genome_mt.rows()[exome_mt.row_key]
        exome_mt = exome_mt.annotate_rows(
            exome_pass=hl.len(exome_mt.filters) == 0,
            genome_pass=hl.len(genome_rows_join.filters) == 0,
            genome_freq=genome_rows_join.freq
        )

        def join_struct(exome_field, join_field, exome_an_field, genome_an_field):
            ac = exome_field.AC[1] + hl.or_else(join_field.AC[1], 0)
            an = exome_an_field + hl.or_else(genome_an_field, 0)
            hom_count = exome_field.homozygote_count[1] + hl.or_else(join_field.homozygote_count[1], 0)
            af = ac / an
            expected_homs = an * (af ** 2) / 2
            bi = hom_count / expected_homs
            return hl.struct(AC=ac, AN=an, AF=af, homozygote_count=hom_count, expected_homs=expected_homs, bi=bi)

        def format_pop_name_hail(pop):
            return (hl.case()
                    .when(pop == 'all', 'gnomad')
                    .when(pop == 'fin', 'gnomad_fin')
                    .default('gnomad_nfe_' + pop))

        genome_entries = genome_mt[exome_mt.row_key, exome_mt.col_key]
        genome_index = genome_mt.index_globals().freq_index_dict.get(format_pop_name_hail(exome_mt.pop))
        genome_an = hl.cond(hl.is_defined(genome_index), exome_mt.genome_freq[genome_index].AN, 0)
        mt = exome_mt.annotate_entries(
            **{ann: join_struct(
                exome_mt[ann], genome_entries[ann],
                exome_mt.freq[exome_mt.freq_index_dict[format_pop_name_hail(exome_mt.pop)]].AN, genome_an)
                for ann in annotations.keys()},
        )
        def agg_struct(ac, an, homozygote_count):
            ac = hl.int32(hl.agg.sum(ac))
            an = hl.int32(hl.agg.sum(an))
            hom_count = hl.int32(hl.agg.sum(homozygote_count))
            af = ac / an
            expected_homs = an * (af ** 2) / 2
            bi = hom_count / expected_homs
            return hl.struct(AC=ac, AN=an, AF=af, homozygote_count=hom_count, expected_homs=expected_homs, bi=bi)

        full_mt = mt.group_cols_by(pop='all').aggregate(
            **{ann: agg_struct(mt[ann].AC, mt[ann].AN, mt[ann].homozygote_count)
            for ann in annotations.keys()}
        )
        mt = mt.union_cols(full_mt)
        mt = mt.annotate_cols(ccr5_result=ccr5_freqs[mt.pop])
        mt = mt.filter_rows(mt.locus.in_autosome())

        def compute_p_against_background(mt, tol = 0.01):
            reported_criteria = ((mt.reported.AF >= mt.ccr5_result.reported.AF[1] - tol) &
                                               (mt.reported.AF <= mt.ccr5_result.reported.AF[1] + tol) &
                                               (mt.reported.AN > 0.8 * mt.ccr5_result.reported.AN))
            recomputed_criteria = ((mt.recomputed.AF >= mt.ccr5_result.recomputed.AF[1] - tol) &
                                             (mt.recomputed.AF <= mt.ccr5_result.recomputed.AF[1] + tol) &
                                             (mt.recomputed.AN > 0.8 * mt.ccr5_result.recomputed.AN))
            recomputed_criteria_pass = ((mt.recomputed.AF >= mt.ccr5_result.recomputed.AF[1] - tol) &
                                    (mt.recomputed.AF <= mt.ccr5_result.recomputed.AF[1] + tol) & 
                                    mt.exome_pass & mt.genome_pass &
                                    (mt.recomputed.AN > 0.8 * mt.ccr5_result.recomputed.AN))

            return mt.annotate_cols(
                # Reported variant selection
                p_reported_callrate=hl.agg.filter(reported_criteria, hl.agg.fraction(mt.reported.bi < mt.ccr5_result.reported_bi)),
                total_reported_callrate=hl.agg.count_where(reported_criteria),
                mean_bi_reported_callrate=hl.agg.filter(reported_criteria, hl.agg.mean(mt.reported.bi)),
                median_bi_reported_callrate=hl.agg.filter(reported_criteria, hl.agg.approx_quantiles(mt.reported.bi, 0.5)),

                # Redone variant selection
                p_recomputed_callrate=hl.agg.filter(recomputed_criteria, hl.agg.fraction(mt.recomputed.bi < mt.ccr5_result.recomputed_bi)),
                total_recomputed_callrate=hl.agg.count_where(recomputed_criteria),
                mean_bi_recomputed_callrate=hl.agg.filter(recomputed_criteria, hl.agg.mean(mt.recomputed.bi)),
                median_bi_recomputed_callrate=hl.agg.filter(recomputed_criteria, hl.agg.approx_quantiles(mt.recomputed.bi, 0.5)),

                # Redone variants selection PASS
                p_recomputed_pass=hl.agg.filter(recomputed_criteria_pass, hl.agg.fraction(mt.recomputed.bi < mt.ccr5_result.recomputed_bi)),
                total_recomputed_pass=hl.agg.count_where(recomputed_criteria_pass),
                mean_bi_recomputed_pass=hl.agg.filter(recomputed_criteria_pass, hl.agg.mean(mt.recomputed.bi)),
                median_bi_recomputed_pass=hl.agg.filter(recomputed_criteria_pass, hl.agg.approx_quantiles(mt.recomputed.bi, 0.5)),

            ).cols()
        compute_p_against_background(mt).write(f'{temp_root}/ccr5_recalc.ht', overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--find_high_ad_hets', help='Find high AD hets from public data', action='store_true')
    parser.add_argument('--find_high_ad_hets_by_sample', help='Find high AD hets from individual level data', action='store_true')
    parser.add_argument('--recompute_ccr5_homozygosity', help='Recompute frequency and HWE for CCR5-del32', action='store_true')
    parser.add_argument('--remove_contaminated', help='Remove contaminated samples from CCR5 computation', action='store_true')
    parser.add_argument('--recompute_background_distribution', help='Recompute frequency and HWE for CCR5-del32', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite all data', action='store_true')
    args = parser.parse_args()

    main(args)
