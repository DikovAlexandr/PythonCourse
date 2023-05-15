from PIL import Image
import os

file_list = os.listdir('maps/')
image_files = [file for file in file_list]
image_files.sort()

tau_values = []
images = []
for filename in image_files:
    im = Image.open(f'maps/{filename}')
    images.append(im)

# Сохраняем gif-изображение
images[0].save('maps_by_exe.gif', save_all=True, append_images=images[1:], duration=500, loop=0)