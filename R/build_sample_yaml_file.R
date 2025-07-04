output_yaml = "/home/s.gregoricchio/DATA_seb/Organoids/ChIP/H3K27ac/sample_config.yaml"
config_table = "/home/s.gregoricchio/DATA_seb/Organoids/ChIP/H3K27ac/sample_config_table.txt"

config_table = data.table::fread(config_table)



write(x = "chip_dict:\n", file = output_yaml)

for (i in 1:nrow(config_table)) {
  write(file = output_yaml, append = T,
        paste0("  ", config_table$Sample[i], ":\n",
               "    control: ", config_table$Control[i],"\n",
               "    broad: ", stringr::str_to_sentence(config_table$Broad[i]), "\n"))
}