#batch running to prepare gvcf as ml ready
project="project-FF5gKb80KQgQ2ZbgKZ7GFPYj"
applet="applet-FP0GzpQ0KQgy60j61ZxbFqgp"
address="1000G_gvcf"
output="$address/mltrain"

dx select $project

for i in $(dx find data --name "*.g.vcf.gz" --path $address --brief); do
    dx describe $i --name
    dx run $applet \
    -iInput_array=$i \
    -iInput_array=file-FK4k2J80b4f074KyGgQ3fFBX \
    -iInput_array=file-FY2z1x00yp626BBZ5jjKZ7vJ \
    -iInput_array=file-ByGVY4j04Y0YJ0yJpj0f8qPG \
    -iInput_array=file-Fx03vFj05jBkV6KGPX2G1BxX \
    -iInput_array=file-Fx03vFj05jBQFfp0FKVk4x03 \
    -iInput_array=file-FX02GZj0xfkq5p8zJPxPbjqZ \
    -iInput_array=file-Fx3vj200KQgXyq6b5QVj4J7Y \
    -icmd_sh=file-Fx813300KQgpQp0v4Qjqfxgz \
    --destination $output --brief -y

done
