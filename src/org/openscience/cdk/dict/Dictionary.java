/* $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 *
 * Copyright (C) 2003  The Chemistry Development Kit (CDK) project
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */
package org.openscience.cdk.dict;

import org.openscience.cdk.tools.LoggingTool;
import java.util.Hashtable;
import org.xml.sax.helpers.*;
import org.xml.sax.*;
import java.io.*;
import java.util.Enumeration;

/**
 * Dictionary with entries.
 *
 * <p>FIXME: this should be replace by a uptodate Dictionary Schema
 * DOM type thing.
 *
 * @author     Egon Willighagen
 * @created    2003-08-23
 * @keyword    dictionary
 */
public class Dictionary {

    private Hashtable entries;
    
    public Dictionary() {
        entries = new Hashtable();
    }
    
    /**
     * Initializes this reader.
     */
    private void initUnmarshall() {
    }
    
    public static Dictionary unmarshal(Reader reader) {
        LoggingTool logger = new LoggingTool("org.openscience.cdk.dict.Dictionary");
        DictionaryHandler handler = new DictionaryHandler();
        XMLReader parser = null;
        try {
            parser = new gnu.xml.aelfred2.XmlReader();
            logger.info("Using Aelfred2 XML parser.");
        } catch (Exception e) {
            logger.error("Could not instantiate Aelfred2 XML reader!");
            logger.debug(e);
        }
        try {
            parser.setFeature("http://xml.org/sax/features/validation", false);
            logger.info("Deactivated validation");
        } catch (SAXException e) {
            logger.warn("Cannot deactivate validation.");
            logger.debug(e);
        }
        parser.setContentHandler(handler);
        Dictionary dict = null;
        try {
            parser.parse(new InputSource(reader));
            dict = handler.getDictionary();
        } catch (IOException e) {
            logger.error("IOException: " + e.toString());
        } catch (SAXException saxe) {
            logger.error("SAXException: " + saxe.getClass().getName());
            logger.debug(saxe);
        }
        return dict;
    }
    
    public void addEntry(Entry entry) {
        entries.put(entry.getID(), entry);
    }
    
    public Entry[] getEntry() {
        int size = entries.size();
        Entry[] entryArray = new Entry[size];
        Enumeration elements = entries.elements();
        int counter = 0;
        while (elements.hasMoreElements() && counter < size) {
            entryArray[counter] = (Entry)elements.nextElement();
            counter++;
        }
        return entryArray;
    }
    
    public boolean hasEntry(String id) {
        return entries.containsKey(id);
    }
}
