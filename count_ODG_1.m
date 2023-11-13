function ODG=count_ODG_1(original_wav,watermarked_wav)
    cd audio
    ODG=PQevalAudio (original_wav, watermarked_wav);
    cd ..
end