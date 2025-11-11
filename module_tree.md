app.R
│
├── left_card (patient_filter_module)
│   └── module: module_left_card.R
│
├── right_card (drug_filter_module)
│   └── module: module_right_card.R
│       └── bnf_table (bnf_table_module)
│
└── middle_card (polyph_module)
    └── module: module_middle_card.R
        │
        ├── module_ce (ce_module) - "Outcome explorer"
        │   └── module: module_ce.R
        │       ├── module_overview (overview_module)
        │       ├── module_prescription_explorer (prescription_explorer_module)
        │       │   ├── module_presc_tto (tto_module)
        │       │   └── module_acute_presc (acute_presc_module)
        │       └── module_ltc_explorer (ltc_explorer_module)
        │
        ├── module_ccm (ccm_module) - "Case-control Matching"
        │   └── module: module_ccm.R
        │
        └── module_cca (cca_module) - "Case-control Analysis"
            └── module: module_cca.R
                ├── module_cca_sociodemographics (sociodemog)
                └── module_cca_prevalence (prevalence)