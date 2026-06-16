import io
import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import clonalyzer


BASE_HEADER = (
    "t_hr,Clone,Rep,is_post_feed,VCD,DCD,Viab_pct,rP_mg_L,Glc_g_L,"
    "Lac_g_L,Gln_mM,Glu_mM,Vol_mL,GFP_mean\n"
)


def test_fedbatch_without_volume_falls_back_to_batch_with_warning():
    csv_text = (
        "t_hr,Clone,Rep,is_post_feed,VCD,DCD,Viab_pct,rP_mg_L,Glc_g_L,Lac_g_L,Gln_mM,Glu_mM\n"
        "0,A,1,FALSE,1.0e6,1.0e4,99,0,5.0,0.2,4.0,0.5\n"
        "24,A,1,TRUE,2.0e6,2.0e4,98,10,4.0,0.8,3.0,1.0\n"
    )
    result = clonalyzer.run_analysis(csv_text, skip_first_row=False, use_volume=True)

    assert result["info"]["requested_mode"] == "fedbatch"
    assert result["info"]["effective_mode"] == "batch"
    assert result["info"]["scenario"] == clonalyzer.SCENARIO_CONST
    assert any("Vol_mL" in warning for warning in result["info"]["mode_warnings"])


def test_fedbatch_without_feed_column_falls_back_to_batch_with_warning():
    csv_text = (
        "t_hr,Clone,Rep,VCD,DCD,Viab_pct,rP_mg_L,Glc_g_L,Lac_g_L,Gln_mM,Glu_mM,Vol_mL\n"
        "0,A,1,1.0e6,1.0e4,99,0,5.0,0.2,4.0,0.5,50\n"
        "24,A,1,2.0e6,2.0e4,98,10,4.0,0.8,3.0,1.0,55\n"
    )
    result = clonalyzer.run_analysis(csv_text, skip_first_row=False, use_volume=True)

    assert result["info"]["effective_mode"] == "batch"
    assert any("is_post_feed" in warning for warning in result["info"]["mode_warnings"])


def test_fedbatch_without_true_postfeed_falls_back_to_batch_with_warning():
    csv_text = (
        BASE_HEADER +
        "0,A,1,FALSE,1.0e6,1.0e4,99,0,5.0,0.2,4.0,0.5,50,100\n"
        "24,A,1,FALSE,2.0e6,2.0e4,98,10,4.0,0.8,3.0,1.0,55,\n"
    )
    result = clonalyzer.run_analysis(csv_text, skip_first_row=False, use_volume=True)

    assert result["info"]["effective_mode"] == "batch"
    assert any("is_post_feed = TRUE" in warning for warning in result["info"]["mode_warnings"])


def test_valid_fedbatch_uses_variable_volume_mode():
    csv_text = (
        BASE_HEADER +
        "0,A,1,FALSE,1.0e6,1.0e4,99,0,5.0,0.2,4.0,0.5,50,100\n"
        "24,A,1,TRUE,2.0e6,2.0e4,98,10,4.5,0.8,3.0,1.0,55,\n"
        "48,A,1,TRUE,3.0e6,3.0e4,97,20,4.0,1.2,2.5,1.3,60,120\n"
    )
    result = clonalyzer.run_analysis(csv_text, skip_first_row=False, use_volume=True)

    assert result["info"]["effective_mode"] == "fedbatch"
    assert result["info"]["scenario"] == clonalyzer.SCENARIO_VAR
    assert result["info"]["mode_warnings"] == []


def test_missing_fluorescence_is_not_forward_filled_or_used_for_correlation():
    csv_text = (
        BASE_HEADER +
        "0,A,1,FALSE,1.0e6,1.0e4,99,0,5.0,0.2,4.0,0.5,50,100\n"
        "24,A,1,FALSE,2.0e6,2.0e4,98,10,4.0,0.8,3.0,1.0,50,\n"
    )
    result = clonalyzer.run_analysis(csv_text, skip_first_row=False, use_volume=False)
    processed = pd.read_csv(io.StringIO(result["processed_csv"]))
    corr = clonalyzer.make_custom_correlation("GFP_mean", "rP_mg_L")

    assert pd.isna(processed.loc[1, "GFP_mean"])
    assert corr is not None
    assert len(corr["traces"]) == 1
    assert corr["traces"][0]["x"] == [100.0]
