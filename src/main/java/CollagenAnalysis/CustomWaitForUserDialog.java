package CollagenAnalysis;

import ij.gui.WaitForUserDialog;

import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

public class CustomWaitForUserDialog  extends WaitForUserDialog {

    private boolean escPressed = false;
    public CustomWaitForUserDialog(String title, String text) {
        super(title, text);
    }
    private void addEscapeKeyListener() {
        this.addKeyListener(new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent e) {
                if (e.getKeyCode() == KeyEvent.VK_ESCAPE) {
                    escPressed = true;
                    dispose(); // Close the dialog when Escape is pressed
                }
            }
        });
    }

    public boolean wasEscPressed() {
        return escPressed;
    }
}
