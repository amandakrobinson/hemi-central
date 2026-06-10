# render images
from psychopy import visual
import os
from glob import glob

os.makedirs('stimuli_rendered', exist_ok=True)

# Create window at a reasonable size for rendering
win = visual.Window([384, 384], units='pix', color=(-.1, -.1, -.1))

stimsize = 128*3
wordsize = 40*3

# Object stimuli
objectstimuli = sorted(glob('stims/*.png'))
print(f"Found {len(objectstimuli)} images:")
for img in objectstimuli:
    print(f"  {img}")

for i, imgpath in enumerate(objectstimuli):
    print(f"Processing image {i}: {imgpath}")
    try:
        stim = visual.ImageStim(win, imgpath, size=stimsize)
        stim.draw()
        win.flip()
        img = win.getMovieFrame(buffer='front')
        filename = imgpath.split('/')[-1].replace('.png', '')
        output_path = f'stimuli_rendered/{filename}.png'
        img.save(output_path)
        print(f"  Saved to {output_path}")
    except Exception as e:
        print(f"  ERROR: {e}")

# Word stimuli
wordstimuli = ['fish','bird','face','boat','tree','tool']
word_idx = 25

for word in wordstimuli:
    variants = [
        (word.lower(), 'lower'),
        (word.upper(), 'upper'),
    ]
    
    for text, suffix in variants:
        stim = visual.TextStim(win, text=text, height=wordsize, color='white')
        stim.draw()
        win.flip()
        img = win.getMovieFrame(buffer='front')
        img.save(f'stimuli_rendered/stim0{word_idx:02d}_{word}_{suffix}.png')
        word_idx += 1

win.close()