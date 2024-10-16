
#!/bin/bash

cd /nexus/posix0/MAGE-flaski/service/projects/data/Bioinformatics/bit_gpt_age/
mkdir -p pmc_papers
cd pmc_papers

#################


wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.incr.2023-08-02.tar.gz
tar -zxvf author_manuscript_txt.incr.2023-08-02.tar.gz PMC010xxxxxx/PMC10388377.txt
rm -rf author_manuscript_txt.incr.2023-08-02.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_other/txt/oa_other_txt.PMC004xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_other_txt.PMC004xxxxxx.baseline.2023-06-17.tar.gz PMC004xxxxxx/PMC4858780.txt PMC004xxxxxx/PMC4250049.txt PMC004xxxxxx/PMC4338562.txt PMC004xxxxxx/PMC4161471.txt
rm -rf oa_other_txt.PMC004xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.incr.2023-07-12.tar.gz
tar -zxvf oa_comm_txt.incr.2023-07-12.tar.gz PMC010xxxxxx/PMC10328083.txt
rm -rf oa_comm_txt.incr.2023-07-12.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.PMC008xxxxxx.baseline.2023-06-16.tar.gz
tar -zxvf author_manuscript_txt.PMC008xxxxxx.baseline.2023-06-16.tar.gz PMC008xxxxxx/PMC8259782.txt PMC008xxxxxx/PMC8460615.txt PMC008xxxxxx/PMC8785908.txt PMC008xxxxxx/PMC8592384.txt
rm -rf author_manuscript_txt.PMC008xxxxxx.baseline.2023-06-16.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.PMC004xxxxxx.baseline.2023-06-16.tar.gz
tar -zxvf author_manuscript_txt.PMC004xxxxxx.baseline.2023-06-16.tar.gz PMC004xxxxxx/PMC4445605.txt PMC004xxxxxx/PMC4338562.txt PMC004xxxxxx/PMC4445606.txt PMC004xxxxxx/PMC4016105.txt PMC004xxxxxx/PMC4154513.txt PMC004xxxxxx/PMC4335172.txt PMC004xxxxxx/PMC4427240.txt PMC004xxxxxx/PMC4773887.txt PMC004xxxxxx/PMC4250049.txt PMC004xxxxxx/PMC4684691.txt PMC004xxxxxx/PMC4123609.txt PMC004xxxxxx/PMC4006270.txt PMC004xxxxxx/PMC4547605.txt PMC004xxxxxx/PMC4789495.txt PMC004xxxxxx/PMC4944841.txt PMC004xxxxxx/PMC4051997.txt PMC004xxxxxx/PMC4335171.txt PMC004xxxxxx/PMC4943223.txt PMC004xxxxxx/PMC4161471.txt PMC004xxxxxx/PMC4757735.txt PMC004xxxxxx/PMC4040421.txt PMC004xxxxxx/PMC4467031.txt PMC004xxxxxx/PMC4346203.txt
rm -rf author_manuscript_txt.PMC004xxxxxx.baseline.2023-06-16.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.incr.2023-06-20.tar.gz
tar -zxvf oa_comm_txt.incr.2023-06-20.tar.gz PMC010xxxxxx/PMC10272782.txt
rm -rf oa_comm_txt.incr.2023-06-20.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_other/txt/oa_other_txt.PMC006xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_other_txt.PMC006xxxxxx.baseline.2023-06-17.tar.gz PMC006xxxxxx/PMC6596419.txt PMC006xxxxxx/PMC6713615.txt PMC006xxxxxx/PMC6090564.txt PMC006xxxxxx/PMC6861129.txt PMC006xxxxxx/PMC6225988.txt
rm -rf oa_other_txt.PMC006xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.incr.2023-08-11.tar.gz
tar -zxvf oa_comm_txt.incr.2023-08-11.tar.gz PMC010xxxxxx/PMC10410055.txt
rm -rf oa_comm_txt.incr.2023-08-11.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.incr.2023-07-17.tar.gz
tar -zxvf author_manuscript_txt.incr.2023-07-17.tar.gz PMC004xxxxxx/PMC4684691.txt
rm -rf author_manuscript_txt.incr.2023-07-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_other/txt/oa_other_txt.PMC001xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_other_txt.PMC001xxxxxx.baseline.2023-06-17.tar.gz PMC001xxxxxx/PMC1994694.txt
rm -rf oa_other_txt.PMC001xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.incr.2023-07-30.tar.gz
tar -zxvf oa_comm_txt.incr.2023-07-30.tar.gz PMC010xxxxxx/PMC10374564.txt
rm -rf oa_comm_txt.incr.2023-07-30.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/txt/oa_noncomm_txt.PMC004xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_noncomm_txt.PMC004xxxxxx.baseline.2023-06-17.tar.gz PMC004xxxxxx/PMC4383974.txt PMC004xxxxxx/PMC4937333.txt PMC004xxxxxx/PMC4261080.txt PMC004xxxxxx/PMC4680133.txt PMC004xxxxxx/PMC4664167.txt PMC004xxxxxx/PMC4215182.txt PMC004xxxxxx/PMC4896057.txt PMC004xxxxxx/PMC4052133.txt PMC004xxxxxx/PMC4431669.txt PMC004xxxxxx/PMC4941146.txt PMC004xxxxxx/PMC4975551.txt PMC004xxxxxx/PMC4332246.txt PMC004xxxxxx/PMC4267620.txt PMC004xxxxxx/PMC4471822.txt PMC004xxxxxx/PMC4138322.txt PMC004xxxxxx/PMC4738383.txt PMC004xxxxxx/PMC4559908.txt PMC004xxxxxx/PMC4792078.txt PMC004xxxxxx/PMC4339303.txt PMC004xxxxxx/PMC4752821.txt PMC004xxxxxx/PMC4033770.txt
rm -rf oa_noncomm_txt.PMC004xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.incr.2023-06-25.tar.gz
tar -zxvf oa_comm_txt.incr.2023-06-25.tar.gz PMC010xxxxxx/PMC10287738.txt
rm -rf oa_comm_txt.incr.2023-06-25.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/txt/oa_noncomm_txt.PMC007xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_noncomm_txt.PMC007xxxxxx.baseline.2023-06-17.tar.gz PMC007xxxxxx/PMC7237863.txt PMC007xxxxxx/PMC7033349.txt PMC007xxxxxx/PMC7455195.txt PMC007xxxxxx/PMC7337487.txt PMC007xxxxxx/PMC7895438.txt PMC007xxxxxx/PMC7036479.txt PMC007xxxxxx/PMC7198979.txt PMC007xxxxxx/PMC7868890.txt PMC007xxxxxx/PMC7361206.txt PMC007xxxxxx/PMC7383232.txt PMC007xxxxxx/PMC7567641.txt PMC007xxxxxx/PMC7560164.txt PMC007xxxxxx/PMC7272126.txt
rm -rf oa_noncomm_txt.PMC007xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.incr.2023-07-13.tar.gz
tar -zxvf author_manuscript_txt.incr.2023-07-13.tar.gz PMC004xxxxxx/PMC4547605.txt
rm -rf author_manuscript_txt.incr.2023-07-13.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/txt/oa_noncomm_txt.PMC005xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_noncomm_txt.PMC005xxxxxx.baseline.2023-06-17.tar.gz PMC005xxxxxx/PMC5223084.txt PMC005xxxxxx/PMC5886115.txt PMC005xxxxxx/PMC5721940.txt PMC005xxxxxx/PMC5584192.txt PMC005xxxxxx/PMC5437728.txt PMC005xxxxxx/PMC5933446.txt PMC005xxxxxx/PMC5628164.txt PMC005xxxxxx/PMC5947121.txt PMC005xxxxxx/PMC5001613.txt PMC005xxxxxx/PMC5389577.txt PMC005xxxxxx/PMC5071710.txt PMC005xxxxxx/PMC5144889.txt PMC005xxxxxx/PMC5499872.txt PMC005xxxxxx/PMC5940314.txt
rm -rf oa_noncomm_txt.PMC005xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.PMC009xxxxxx.baseline.2023-06-16.tar.gz
tar -zxvf author_manuscript_txt.PMC009xxxxxx.baseline.2023-06-16.tar.gz PMC009xxxxxx/PMC9999291.txt PMC009xxxxxx/PMC9027952.txt PMC009xxxxxx/PMC9443067.txt
rm -rf author_manuscript_txt.PMC009xxxxxx.baseline.2023-06-16.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/txt/oa_noncomm_txt.PMC009xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_noncomm_txt.PMC009xxxxxx.baseline.2023-06-17.tar.gz PMC009xxxxxx/PMC9676392.txt PMC009xxxxxx/PMC9044173.txt PMC009xxxxxx/PMC9873946.txt PMC009xxxxxx/PMC9016346.txt PMC009xxxxxx/PMC9499549.txt PMC009xxxxxx/PMC9027952.txt PMC009xxxxxx/PMC9727813.txt
rm -rf oa_noncomm_txt.PMC009xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/txt/oa_noncomm_txt.PMC003xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_noncomm_txt.PMC003xxxxxx.baseline.2023-06-17.tar.gz PMC003xxxxxx/PMC3330212.txt PMC003xxxxxx/PMC3633371.txt PMC003xxxxxx/PMC3973307.txt PMC003xxxxxx/PMC3018165.txt PMC003xxxxxx/PMC3781460.txt PMC003xxxxxx/PMC3749354.txt PMC003xxxxxx/PMC3875664.txt PMC003xxxxxx/PMC3988997.txt PMC003xxxxxx/PMC3902947.txt PMC003xxxxxx/PMC3105550.txt PMC003xxxxxx/PMC3019561.txt PMC003xxxxxx/PMC3998800.txt PMC003xxxxxx/PMC3543749.txt PMC003xxxxxx/PMC3328357.txt PMC003xxxxxx/PMC3159968.txt PMC003xxxxxx/PMC3439893.txt
rm -rf oa_noncomm_txt.PMC003xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/txt/oa_noncomm_txt.PMC010xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_noncomm_txt.PMC010xxxxxx.baseline.2023-06-17.tar.gz PMC010xxxxxx/PMC10173804.txt
rm -rf oa_noncomm_txt.PMC010xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.PMC005xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_comm_txt.PMC005xxxxxx.baseline.2023-06-17.tar.gz PMC005xxxxxx/PMC5582349.txt PMC005xxxxxx/PMC5691074.txt PMC005xxxxxx/PMC5580880.txt PMC005xxxxxx/PMC5457527.txt PMC005xxxxxx/PMC5416710.txt PMC005xxxxxx/PMC5760849.txt PMC005xxxxxx/PMC5032696.txt PMC005xxxxxx/PMC5656751.txt PMC005xxxxxx/PMC5370449.txt PMC005xxxxxx/PMC5577152.txt PMC005xxxxxx/PMC5059774.txt PMC005xxxxxx/PMC5983821.txt PMC005xxxxxx/PMC5333540.txt PMC005xxxxxx/PMC5536874.txt PMC005xxxxxx/PMC5242300.txt PMC005xxxxxx/PMC5095294.txt PMC005xxxxxx/PMC5395972.txt PMC005xxxxxx/PMC5837557.txt PMC005xxxxxx/PMC5595691.txt PMC005xxxxxx/PMC5426901.txt PMC005xxxxxx/PMC5772647.txt PMC005xxxxxx/PMC5933513.txt PMC005xxxxxx/PMC5087735.txt PMC005xxxxxx/PMC5135034.txt PMC005xxxxxx/PMC5627050.txt PMC005xxxxxx/PMC5786504.txt PMC005xxxxxx/PMC5434020.txt PMC005xxxxxx/PMC5582868.txt PMC005xxxxxx/PMC5783719.txt PMC005xxxxxx/PMC5425116.txt PMC005xxxxxx/PMC5155225.txt PMC005xxxxxx/PMC5839405.txt PMC005xxxxxx/PMC5992170.txt PMC005xxxxxx/PMC5597081.txt PMC005xxxxxx/PMC5552275.txt PMC005xxxxxx/PMC5820258.txt PMC005xxxxxx/PMC5055666.txt PMC005xxxxxx/PMC5610247.txt PMC005xxxxxx/PMC5865154.txt PMC005xxxxxx/PMC5402976.txt PMC005xxxxxx/PMC5344681.txt PMC005xxxxxx/PMC5507281.txt PMC005xxxxxx/PMC5355364.txt PMC005xxxxxx/PMC5834439.txt PMC005xxxxxx/PMC5835776.txt PMC005xxxxxx/PMC5665913.txt PMC005xxxxxx/PMC5026704.txt PMC005xxxxxx/PMC5240381.txt PMC005xxxxxx/PMC5310665.txt PMC005xxxxxx/PMC5013015.txt PMC005xxxxxx/PMC5735189.txt PMC005xxxxxx/PMC5988561.txt PMC005xxxxxx/PMC5032245.txt PMC005xxxxxx/PMC5703644.txt PMC005xxxxxx/PMC5749727.txt PMC005xxxxxx/PMC5333801.txt PMC005xxxxxx/PMC5039181.txt PMC005xxxxxx/PMC5566455.txt PMC005xxxxxx/PMC5615923.txt
rm -rf oa_comm_txt.PMC005xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_other/txt/oa_other_txt.PMC008xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_other_txt.PMC008xxxxxx.baseline.2023-06-17.tar.gz PMC008xxxxxx/PMC8102037.txt
rm -rf oa_other_txt.PMC008xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_other/txt/oa_other_txt.PMC005xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_other_txt.PMC005xxxxxx.baseline.2023-06-17.tar.gz PMC005xxxxxx/PMC5113755.txt PMC005xxxxxx/PMC5088783.txt PMC005xxxxxx/PMC5521867.txt
rm -rf oa_other_txt.PMC005xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.PMC002xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_comm_txt.PMC002xxxxxx.baseline.2023-06-17.tar.gz PMC002xxxxxx/PMC2683749.txt PMC002xxxxxx/PMC2932684.txt PMC002xxxxxx/PMC2942829.txt PMC002xxxxxx/PMC2698149.txt PMC002xxxxxx/PMC2824086.txt PMC002xxxxxx/PMC2258177.txt PMC002xxxxxx/PMC2214808.txt PMC002xxxxxx/PMC2829060.txt PMC002xxxxxx/PMC2098796.txt PMC002xxxxxx/PMC2806523.txt PMC002xxxxxx/PMC2841611.txt PMC002xxxxxx/PMC2614481.txt
rm -rf oa_comm_txt.PMC002xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.incr.2023-06-21.tar.gz
tar -zxvf author_manuscript_txt.incr.2023-06-21.tar.gz PMC009xxxxxx/PMC9027952.txt PMC003xxxxxx/PMC3840926.txt
rm -rf author_manuscript_txt.incr.2023-06-21.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.incr.2023-08-02.tar.gz
tar -zxvf oa_comm_txt.incr.2023-08-02.tar.gz PMC010xxxxxx/PMC10388377.txt
rm -rf oa_comm_txt.incr.2023-08-02.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.incr.2023-06-23.tar.gz
tar -zxvf author_manuscript_txt.incr.2023-06-23.tar.gz PMC003xxxxxx/PMC3909615.txt
rm -rf author_manuscript_txt.incr.2023-06-23.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.PMC006xxxxxx.baseline.2023-06-16.tar.gz
tar -zxvf author_manuscript_txt.PMC006xxxxxx.baseline.2023-06-16.tar.gz PMC006xxxxxx/PMC6696993.txt PMC006xxxxxx/PMC6284793.txt PMC006xxxxxx/PMC6596419.txt PMC006xxxxxx/PMC6713615.txt PMC006xxxxxx/PMC6057610.txt PMC006xxxxxx/PMC6090564.txt PMC006xxxxxx/PMC6861129.txt PMC006xxxxxx/PMC6881776.txt PMC006xxxxxx/PMC6245547.txt PMC006xxxxxx/PMC6225988.txt PMC006xxxxxx/PMC6173523.txt PMC006xxxxxx/PMC6899219.txt PMC006xxxxxx/PMC6231542.txt PMC006xxxxxx/PMC6942693.txt PMC006xxxxxx/PMC6814438.txt PMC006xxxxxx/PMC6541012.txt
rm -rf author_manuscript_txt.PMC006xxxxxx.baseline.2023-06-16.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_other/txt/oa_other_txt.PMC007xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_other_txt.PMC007xxxxxx.baseline.2023-06-17.tar.gz PMC007xxxxxx/PMC7834195.txt PMC007xxxxxx/PMC7105183.txt PMC007xxxxxx/PMC7577179.txt
rm -rf oa_other_txt.PMC007xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.PMC007xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_comm_txt.PMC007xxxxxx.baseline.2023-06-17.tar.gz PMC007xxxxxx/PMC7182428.txt PMC007xxxxxx/PMC7971964.txt PMC007xxxxxx/PMC7797037.txt PMC007xxxxxx/PMC7294789.txt PMC007xxxxxx/PMC7932915.txt PMC007xxxxxx/PMC7058826.txt PMC007xxxxxx/PMC7910479.txt PMC007xxxxxx/PMC7769605.txt PMC007xxxxxx/PMC7803558.txt PMC007xxxxxx/PMC7936346.txt PMC007xxxxxx/PMC7672635.txt PMC007xxxxxx/PMC7266280.txt PMC007xxxxxx/PMC7366647.txt PMC007xxxxxx/PMC7649373.txt PMC007xxxxxx/PMC7138244.txt PMC007xxxxxx/PMC7412830.txt PMC007xxxxxx/PMC7490873.txt PMC007xxxxxx/PMC7491019.txt PMC007xxxxxx/PMC7949116.txt PMC007xxxxxx/PMC7960713.txt PMC007xxxxxx/PMC7602196.txt PMC007xxxxxx/PMC7435156.txt PMC007xxxxxx/PMC7771071.txt PMC007xxxxxx/PMC7101016.txt PMC007xxxxxx/PMC7519657.txt PMC007xxxxxx/PMC7287080.txt PMC007xxxxxx/PMC7317364.txt PMC007xxxxxx/PMC7859412.txt PMC007xxxxxx/PMC7738165.txt PMC007xxxxxx/PMC7320053.txt PMC007xxxxxx/PMC7030929.txt PMC007xxxxxx/PMC7860130.txt PMC007xxxxxx/PMC7332298.txt PMC007xxxxxx/PMC7576295.txt PMC007xxxxxx/PMC7901890.txt PMC007xxxxxx/PMC7379855.txt PMC007xxxxxx/PMC7287348.txt PMC007xxxxxx/PMC7000685.txt PMC007xxxxxx/PMC7338799.txt PMC007xxxxxx/PMC7462614.txt PMC007xxxxxx/PMC7253064.txt PMC007xxxxxx/PMC7318944.txt PMC007xxxxxx/PMC7424762.txt PMC007xxxxxx/PMC7504350.txt PMC007xxxxxx/PMC7649784.txt PMC007xxxxxx/PMC7923679.txt PMC007xxxxxx/PMC7895802.txt
rm -rf oa_comm_txt.PMC007xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.PMC002xxxxxx.baseline.2023-06-16.tar.gz
tar -zxvf author_manuscript_txt.PMC002xxxxxx.baseline.2023-06-16.tar.gz PMC002xxxxxx/PMC2757405.txt PMC002xxxxxx/PMC2659505.txt PMC002xxxxxx/PMC2939271.txt PMC002xxxxxx/PMC2989461.txt PMC002xxxxxx/PMC2902793.txt PMC002xxxxxx/PMC2712670.txt PMC002xxxxxx/PMC2744646.txt
rm -rf author_manuscript_txt.PMC002xxxxxx.baseline.2023-06-16.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/txt/oa_noncomm_txt.PMC008xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_noncomm_txt.PMC008xxxxxx.baseline.2023-06-17.tar.gz PMC008xxxxxx/PMC8933812.txt PMC008xxxxxx/PMC8986112.txt PMC008xxxxxx/PMC8108515.txt PMC008xxxxxx/PMC8164091.txt PMC008xxxxxx/PMC8906734.txt PMC008xxxxxx/PMC8020789.txt PMC008xxxxxx/PMC8097689.txt PMC008xxxxxx/PMC8589318.txt PMC008xxxxxx/PMC8382303.txt PMC008xxxxxx/PMC8168506.txt PMC008xxxxxx/PMC8480666.txt PMC008xxxxxx/PMC8684738.txt PMC008xxxxxx/PMC8462889.txt
rm -rf oa_noncomm_txt.PMC008xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.PMC006xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_comm_txt.PMC006xxxxxx.baseline.2023-06-17.tar.gz PMC006xxxxxx/PMC6677768.txt PMC006xxxxxx/PMC6051193.txt PMC006xxxxxx/PMC6553695.txt PMC006xxxxxx/PMC6702196.txt PMC006xxxxxx/PMC6197506.txt PMC006xxxxxx/PMC6690953.txt PMC006xxxxxx/PMC6408443.txt PMC006xxxxxx/PMC6657895.txt PMC006xxxxxx/PMC6369905.txt PMC006xxxxxx/PMC6609229.txt PMC006xxxxxx/PMC6717989.txt PMC006xxxxxx/PMC6863466.txt PMC006xxxxxx/PMC6203810.txt PMC006xxxxxx/PMC6323124.txt PMC006xxxxxx/PMC6800352.txt PMC006xxxxxx/PMC6469880.txt PMC006xxxxxx/PMC6693550.txt PMC006xxxxxx/PMC6811592.txt PMC006xxxxxx/PMC6909895.txt PMC006xxxxxx/PMC6363466.txt PMC006xxxxxx/PMC6859486.txt PMC006xxxxxx/PMC6172693.txt PMC006xxxxxx/PMC6385858.txt PMC006xxxxxx/PMC6362144.txt PMC006xxxxxx/PMC6599210.txt PMC006xxxxxx/PMC6558167.txt PMC006xxxxxx/PMC6156541.txt PMC006xxxxxx/PMC6811558.txt PMC006xxxxxx/PMC6032845.txt PMC006xxxxxx/PMC6696993.txt PMC006xxxxxx/PMC6127302.txt PMC006xxxxxx/PMC6343102.txt PMC006xxxxxx/PMC6893025.txt PMC006xxxxxx/PMC6096558.txt PMC006xxxxxx/PMC6411718.txt PMC006xxxxxx/PMC6842169.txt PMC006xxxxxx/PMC6914421.txt PMC006xxxxxx/PMC6173523.txt PMC006xxxxxx/PMC6690939.txt PMC006xxxxxx/PMC6694136.txt PMC006xxxxxx/PMC6342327.txt PMC006xxxxxx/PMC6381146.txt PMC006xxxxxx/PMC6448944.txt PMC006xxxxxx/PMC6328943.txt PMC006xxxxxx/PMC6789098.txt PMC006xxxxxx/PMC6051225.txt PMC006xxxxxx/PMC6348486.txt PMC006xxxxxx/PMC6238418.txt PMC006xxxxxx/PMC6514073.txt PMC006xxxxxx/PMC6618114.txt PMC006xxxxxx/PMC6916080.txt PMC006xxxxxx/PMC6974722.txt PMC006xxxxxx/PMC6281273.txt PMC006xxxxxx/PMC6100731.txt PMC006xxxxxx/PMC6755124.txt
rm -rf oa_comm_txt.PMC006xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.PMC008xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_comm_txt.PMC008xxxxxx.baseline.2023-06-17.tar.gz PMC008xxxxxx/PMC8789308.txt PMC008xxxxxx/PMC8782725.txt PMC008xxxxxx/PMC8920897.txt PMC008xxxxxx/PMC8190027.txt PMC008xxxxxx/PMC8553732.txt PMC008xxxxxx/PMC8005083.txt PMC008xxxxxx/PMC8209036.txt PMC008xxxxxx/PMC8654302.txt PMC008xxxxxx/PMC8120713.txt PMC008xxxxxx/PMC8260797.txt PMC008xxxxxx/PMC8217552.txt PMC008xxxxxx/PMC8565073.txt PMC008xxxxxx/PMC8561302.txt PMC008xxxxxx/PMC8881936.txt PMC008xxxxxx/PMC8315803.txt PMC008xxxxxx/PMC8098112.txt PMC008xxxxxx/PMC8024730.txt PMC008xxxxxx/PMC8632672.txt PMC008xxxxxx/PMC8890471.txt PMC008xxxxxx/PMC8987268.txt PMC008xxxxxx/PMC8007214.txt PMC008xxxxxx/PMC8172246.txt PMC008xxxxxx/PMC8728608.txt PMC008xxxxxx/PMC8880994.txt PMC008xxxxxx/PMC8481228.txt PMC008xxxxxx/PMC8421416.txt PMC008xxxxxx/PMC8355499.txt PMC008xxxxxx/PMC8133780.txt PMC008xxxxxx/PMC8412381.txt PMC008xxxxxx/PMC8618896.txt PMC008xxxxxx/PMC8350892.txt PMC008xxxxxx/PMC8548325.txt PMC008xxxxxx/PMC8950555.txt PMC008xxxxxx/PMC8611647.txt PMC008xxxxxx/PMC8185549.txt PMC008xxxxxx/PMC8190293.txt PMC008xxxxxx/PMC8490984.txt PMC008xxxxxx/PMC8520709.txt PMC008xxxxxx/PMC8182411.txt PMC008xxxxxx/PMC8067846.txt PMC008xxxxxx/PMC8460253.txt PMC008xxxxxx/PMC8144018.txt PMC008xxxxxx/PMC8093921.txt PMC008xxxxxx/PMC8719881.txt PMC008xxxxxx/PMC8703237.txt PMC008xxxxxx/PMC8060030.txt PMC008xxxxxx/PMC8408345.txt PMC008xxxxxx/PMC8041777.txt PMC008xxxxxx/PMC8138124.txt PMC008xxxxxx/PMC8465008.txt PMC008xxxxxx/PMC8424529.txt PMC008xxxxxx/PMC8186904.txt PMC008xxxxxx/PMC8050774.txt PMC008xxxxxx/PMC8302645.txt PMC008xxxxxx/PMC8353256.txt PMC008xxxxxx/PMC8789301.txt PMC008xxxxxx/PMC8445024.txt PMC008xxxxxx/PMC8114127.txt PMC008xxxxxx/PMC8970586.txt PMC008xxxxxx/PMC8353526.txt PMC008xxxxxx/PMC8387425.txt PMC008xxxxxx/PMC8163813.txt PMC008xxxxxx/PMC8163834.txt PMC008xxxxxx/PMC8818597.txt
rm -rf oa_comm_txt.PMC008xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.incr.2023-08-04.tar.gz
tar -zxvf author_manuscript_txt.incr.2023-08-04.tar.gz PMC010xxxxxx/PMC10388377.txt
rm -rf author_manuscript_txt.incr.2023-08-04.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/txt/oa_noncomm_txt.PMC002xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_noncomm_txt.PMC002xxxxxx.baseline.2023-06-17.tar.gz PMC002xxxxxx/PMC2063881.txt PMC002xxxxxx/PMC2265410.txt PMC002xxxxxx/PMC2977078.txt PMC002xxxxxx/PMC2775896.txt PMC002xxxxxx/PMC2665216.txt PMC002xxxxxx/PMC2566888.txt PMC002xxxxxx/PMC2119813.txt PMC002xxxxxx/PMC2654118.txt PMC002xxxxxx/PMC2120990.txt PMC002xxxxxx/PMC2367725.txt PMC002xxxxxx/PMC2965244.txt PMC002xxxxxx/PMC2367695.txt PMC002xxxxxx/PMC2806285.txt PMC002xxxxxx/PMC2173682.txt
rm -rf oa_noncomm_txt.PMC002xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.incr.2023-08-03.tar.gz
tar -zxvf oa_comm_txt.incr.2023-08-03.tar.gz PMC010xxxxxx/PMC10391628.txt PMC010xxxxxx/PMC10393417.txt
rm -rf oa_comm_txt.incr.2023-08-03.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.incr.2023-07-20.tar.gz
tar -zxvf oa_comm_txt.incr.2023-07-20.tar.gz PMC010xxxxxx/PMC10352420.txt
rm -rf oa_comm_txt.incr.2023-07-20.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.PMC010xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_comm_txt.PMC010xxxxxx.baseline.2023-06-17.tar.gz PMC010xxxxxx/PMC10148235.txt PMC010xxxxxx/PMC10014586.txt PMC010xxxxxx/PMC10262077.txt PMC010xxxxxx/PMC10154223.txt PMC010xxxxxx/PMC10154239.txt PMC010xxxxxx/PMC10132977.txt PMC010xxxxxx/PMC10205078.txt PMC010xxxxxx/PMC10097726.txt PMC010xxxxxx/PMC10169802.txt PMC010xxxxxx/PMC10154236.txt PMC010xxxxxx/PMC10154229.txt
rm -rf oa_comm_txt.PMC010xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/txt/oa_noncomm_txt.PMC006xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_noncomm_txt.PMC006xxxxxx.baseline.2023-06-17.tar.gz PMC006xxxxxx/PMC6363443.txt PMC006xxxxxx/PMC6315296.txt PMC006xxxxxx/PMC6761770.txt PMC006xxxxxx/PMC6486631.txt PMC006xxxxxx/PMC6549021.txt PMC006xxxxxx/PMC6061787.txt PMC006xxxxxx/PMC6697794.txt PMC006xxxxxx/PMC6900737.txt PMC006xxxxxx/PMC6447380.txt PMC006xxxxxx/PMC6187118.txt PMC006xxxxxx/PMC6745493.txt PMC006xxxxxx/PMC6891098.txt PMC006xxxxxx/PMC6660734.txt
rm -rf oa_noncomm_txt.PMC006xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.incr.2023-07-08.tar.gz
tar -zxvf oa_comm_txt.incr.2023-07-08.tar.gz PMC010xxxxxx/PMC10322686.txt
rm -rf oa_comm_txt.incr.2023-07-08.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.incr.2023-06-23.tar.gz
tar -zxvf oa_comm_txt.incr.2023-06-23.tar.gz PMC010xxxxxx/PMC10281999.txt
rm -rf oa_comm_txt.incr.2023-06-23.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.PMC005xxxxxx.baseline.2023-06-16.tar.gz
tar -zxvf author_manuscript_txt.PMC005xxxxxx.baseline.2023-06-16.tar.gz PMC005xxxxxx/PMC5113755.txt PMC005xxxxxx/PMC5788019.txt PMC005xxxxxx/PMC5415086.txt PMC005xxxxxx/PMC5293008.txt PMC005xxxxxx/PMC5935120.txt PMC005xxxxxx/PMC5507281.txt PMC005xxxxxx/PMC5007171.txt PMC005xxxxxx/PMC5088783.txt
rm -rf author_manuscript_txt.PMC005xxxxxx.baseline.2023-06-16.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.incr.2023-07-04.tar.gz
tar -zxvf oa_comm_txt.incr.2023-07-04.tar.gz PMC010xxxxxx/PMC10272782.txt
rm -rf oa_comm_txt.incr.2023-07-04.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.PMC007xxxxxx.baseline.2023-06-16.tar.gz
tar -zxvf author_manuscript_txt.PMC007xxxxxx.baseline.2023-06-16.tar.gz PMC007xxxxxx/PMC7193988.txt PMC007xxxxxx/PMC7610725.txt PMC007xxxxxx/PMC7603428.txt PMC007xxxxxx/PMC7350392.txt PMC007xxxxxx/PMC7612338.txt PMC007xxxxxx/PMC7610778.txt PMC007xxxxxx/PMC7577179.txt
rm -rf author_manuscript_txt.PMC007xxxxxx.baseline.2023-06-16.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.PMC004xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_comm_txt.PMC004xxxxxx.baseline.2023-06-17.tar.gz PMC004xxxxxx/PMC4056462.txt PMC004xxxxxx/PMC4364945.txt PMC004xxxxxx/PMC4004538.txt PMC004xxxxxx/PMC4531065.txt PMC004xxxxxx/PMC4488357.txt PMC004xxxxxx/PMC4631106.txt PMC004xxxxxx/PMC4992892.txt PMC004xxxxxx/PMC4642670.txt PMC004xxxxxx/PMC4915168.txt PMC004xxxxxx/PMC4270288.txt PMC004xxxxxx/PMC4599244.txt PMC004xxxxxx/PMC4558960.txt PMC004xxxxxx/PMC4767367.txt PMC004xxxxxx/PMC4518038.txt PMC004xxxxxx/PMC4115666.txt PMC004xxxxxx/PMC4191934.txt PMC004xxxxxx/PMC4366498.txt PMC004xxxxxx/PMC4570658.txt PMC004xxxxxx/PMC4141908.txt PMC004xxxxxx/PMC4488845.txt PMC004xxxxxx/PMC4773887.txt PMC004xxxxxx/PMC4557304.txt PMC004xxxxxx/PMC4642562.txt PMC004xxxxxx/PMC4814583.txt PMC004xxxxxx/PMC4360128.txt PMC004xxxxxx/PMC4866704.txt PMC004xxxxxx/PMC4287421.txt PMC004xxxxxx/PMC4518474.txt PMC004xxxxxx/PMC4326868.txt PMC004xxxxxx/PMC4422141.txt PMC004xxxxxx/PMC4770150.txt PMC004xxxxxx/PMC4569844.txt PMC004xxxxxx/PMC4996512.txt PMC004xxxxxx/PMC4382338.txt PMC004xxxxxx/PMC4804169.txt PMC004xxxxxx/PMC4754342.txt PMC004xxxxxx/PMC4850359.txt PMC004xxxxxx/PMC4778021.txt PMC004xxxxxx/PMC4680886.txt PMC004xxxxxx/PMC4991934.txt PMC004xxxxxx/PMC4856093.txt PMC004xxxxxx/PMC4766491.txt PMC004xxxxxx/PMC4169242.txt PMC004xxxxxx/PMC4103672.txt PMC004xxxxxx/PMC4350321.txt PMC004xxxxxx/PMC4653391.txt PMC004xxxxxx/PMC4717275.txt PMC004xxxxxx/PMC4747186.txt PMC004xxxxxx/PMC4222368.txt PMC004xxxxxx/PMC4469414.txt PMC004xxxxxx/PMC4439058.txt PMC004xxxxxx/PMC4515122.txt PMC004xxxxxx/PMC4527767.txt PMC004xxxxxx/PMC4805549.txt PMC004xxxxxx/PMC4546441.txt PMC004xxxxxx/PMC4761717.txt PMC004xxxxxx/PMC4922056.txt PMC004xxxxxx/PMC4740377.txt PMC004xxxxxx/PMC4482076.txt PMC004xxxxxx/PMC4868951.txt PMC004xxxxxx/PMC4018447.txt PMC004xxxxxx/PMC4959029.txt PMC004xxxxxx/PMC4965761.txt PMC004xxxxxx/PMC4770714.txt
rm -rf oa_comm_txt.PMC004xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.PMC009xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_comm_txt.PMC009xxxxxx.baseline.2023-06-17.tar.gz PMC009xxxxxx/PMC9647175.txt PMC009xxxxxx/PMC9833670.txt PMC009xxxxxx/PMC9577955.txt PMC009xxxxxx/PMC9534930.txt PMC009xxxxxx/PMC9009120.txt PMC009xxxxxx/PMC9119528.txt PMC009xxxxxx/PMC9261369.txt PMC009xxxxxx/PMC9135405.txt PMC009xxxxxx/PMC9810596.txt PMC009xxxxxx/PMC9859871.txt PMC009xxxxxx/PMC9984387.txt PMC009xxxxxx/PMC9510030.txt PMC009xxxxxx/PMC9604114.txt PMC009xxxxxx/PMC9206299.txt PMC009xxxxxx/PMC9109056.txt PMC009xxxxxx/PMC9713713.txt PMC009xxxxxx/PMC9881170.txt PMC009xxxxxx/PMC9946356.txt PMC009xxxxxx/PMC9741504.txt PMC009xxxxxx/PMC9410905.txt PMC009xxxxxx/PMC9459415.txt PMC009xxxxxx/PMC9481461.txt PMC009xxxxxx/PMC9841525.txt PMC009xxxxxx/PMC9200640.txt PMC009xxxxxx/PMC9569565.txt PMC009xxxxxx/PMC9600423.txt PMC009xxxxxx/PMC9894012.txt PMC009xxxxxx/PMC9928583.txt PMC009xxxxxx/PMC9741792.txt PMC009xxxxxx/PMC9014876.txt PMC009xxxxxx/PMC9261424.txt PMC009xxxxxx/PMC9909936.txt PMC009xxxxxx/PMC9853457.txt PMC009xxxxxx/PMC9348635.txt PMC009xxxxxx/PMC9344846.txt
rm -rf oa_comm_txt.PMC009xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.PMC003xxxxxx.baseline.2023-06-16.tar.gz
tar -zxvf author_manuscript_txt.PMC003xxxxxx.baseline.2023-06-16.tar.gz PMC003xxxxxx/PMC3836174.txt PMC003xxxxxx/PMC3725963.txt PMC003xxxxxx/PMC3428868.txt PMC003xxxxxx/PMC3909774.txt PMC003xxxxxx/PMC3901671.txt PMC003xxxxxx/PMC3492742.txt PMC003xxxxxx/PMC3918218.txt PMC003xxxxxx/PMC3529641.txt PMC003xxxxxx/PMC3932786.txt PMC003xxxxxx/PMC3909771.txt PMC003xxxxxx/PMC3820420.txt PMC003xxxxxx/PMC3840926.txt PMC003xxxxxx/PMC3624904.txt PMC003xxxxxx/PMC3877687.txt PMC003xxxxxx/PMC3593982.txt PMC003xxxxxx/PMC3924769.txt PMC003xxxxxx/PMC3188402.txt PMC003xxxxxx/PMC3598589.txt PMC003xxxxxx/PMC3909615.txt PMC003xxxxxx/PMC3767426.txt
rm -rf author_manuscript_txt.PMC003xxxxxx.baseline.2023-06-16.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_other/txt/oa_other_txt.PMC003xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_other_txt.PMC003xxxxxx.baseline.2023-06-17.tar.gz PMC003xxxxxx/PMC3189705.txt PMC003xxxxxx/PMC3901671.txt PMC003xxxxxx/PMC3877687.txt PMC003xxxxxx/PMC3492742.txt PMC003xxxxxx/PMC3820420.txt PMC003xxxxxx/PMC3188402.txt PMC003xxxxxx/PMC3063763.txt
rm -rf oa_other_txt.PMC003xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/oa_comm_txt.PMC003xxxxxx.baseline.2023-06-17.tar.gz
tar -zxvf oa_comm_txt.PMC003xxxxxx.baseline.2023-06-17.tar.gz PMC003xxxxxx/PMC3001312.txt PMC003xxxxxx/PMC3395685.txt PMC003xxxxxx/PMC3160332.txt PMC003xxxxxx/PMC3536695.txt PMC003xxxxxx/PMC3444450.txt PMC003xxxxxx/PMC3493444.txt PMC003xxxxxx/PMC3668712.txt PMC003xxxxxx/PMC3967928.txt PMC003xxxxxx/PMC3283687.txt PMC003xxxxxx/PMC3084200.txt PMC003xxxxxx/PMC3567175.txt PMC003xxxxxx/PMC3323522.txt PMC003xxxxxx/PMC3026792.txt PMC003xxxxxx/PMC3323525.txt PMC003xxxxxx/PMC3193374.txt PMC003xxxxxx/PMC3574159.txt PMC003xxxxxx/PMC3817404.txt PMC003xxxxxx/PMC3675642.txt PMC003xxxxxx/PMC3723528.txt PMC003xxxxxx/PMC3120859.txt PMC003xxxxxx/PMC3916286.txt PMC003xxxxxx/PMC3854790.txt PMC003xxxxxx/PMC3176895.txt PMC003xxxxxx/PMC3192828.txt PMC003xxxxxx/PMC3610751.txt PMC003xxxxxx/PMC3751413.txt PMC003xxxxxx/PMC3990160.txt PMC003xxxxxx/PMC3013413.txt PMC003xxxxxx/PMC3087113.txt PMC003xxxxxx/PMC3256132.txt PMC003xxxxxx/PMC3473082.txt PMC003xxxxxx/PMC3519537.txt PMC003xxxxxx/PMC3192837.txt PMC003xxxxxx/PMC3121960.txt PMC003xxxxxx/PMC3972082.txt PMC003xxxxxx/PMC3053371.txt PMC003xxxxxx/PMC3769260.txt PMC003xxxxxx/PMC3998882.txt PMC003xxxxxx/PMC3814316.txt PMC003xxxxxx/PMC3633240.txt PMC003xxxxxx/PMC3969275.txt PMC003xxxxxx/PMC3193372.txt
rm -rf oa_comm_txt.PMC003xxxxxx.baseline.2023-06-17.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/txt/oa_noncomm_txt.incr.2023-07-08.tar.gz
tar -zxvf oa_noncomm_txt.incr.2023-07-08.tar.gz PMC009xxxxxx/PMC9816818.txt
rm -rf oa_noncomm_txt.incr.2023-07-08.tar.gz

#################

wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/author_manuscript_txt.incr.2023-07-27.tar.gz
tar -zxvf author_manuscript_txt.incr.2023-07-27.tar.gz PMC005xxxxxx/PMC5088783.txt
rm -rf author_manuscript_txt.incr.2023-07-27.tar.gz

#################

