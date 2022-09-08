import os

config_dict = {}

# use tex in mpl?
usetex_value = os.getenv("RFE_MPL_text_usetex") not in [ "False", "false", "0", None ]
config_dict.update({"usetex":usetex_value})
