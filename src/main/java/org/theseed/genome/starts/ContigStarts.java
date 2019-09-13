/**
 *
 */
package org.theseed.genome.starts;

import java.util.HashMap;
import java.util.Map;

import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;

/**
 * This class returns all the roles on the specified strand of the contig sorted
 * by start location.  For locations with multiple roles, they are separated by
 * a slash.
 *
 * @author Bruce Parrello
 *
 */
public abstract class ContigStarts {

	// FIELDS
	/** map of locations to role strings */
	Map<Integer,String> locToFunction;

	/**
	 * Create an empty location/role map object.
	 */
	protected ContigStarts() {
		this.locToFunction = new HashMap<Integer, String>();
	}

	/**
	 * @return 	the role string for the specified location, or NULL if there is no
	 * 			feature starting at that location
	 *
	 * @param location	contig location whose roles are desired
	 */
	public String getRoles(int location) {
		return this.locToFunction.get(location);
	}

	/**
	 * Store the role string for a location.
	 *
	 * @param location	relevant contig position
	 * @param roles		string of roles (comma-delimited)
	 */
	protected void setStart(int location, String roles) {
		this.locToFunction.put(location, roles);
	}

	/**
	 * Add a feature to the contig start map.
	 *
	 * @param feat	feature whose start is to be recorded as a real start
	 */
	protected void addFeature(Feature feat) {
		int location = feat.getLocation().getBegin();
		String[] roles = Feature.rolesOfFunction(feat.getFunction());
		// If there is already a start at this location, we simply append
		// more roles.  Otherwise, we start with an empty string.
		StringBuilder roleString = new StringBuilder(20);
		String oldRoles = this.getRoles(location);
		String delim = "";
		if (oldRoles != null) {
			roleString.append(oldRoles);
			delim = " / ";
		}
		// Add the roles found.
		for (String role : roles) {
			String newRole = this.roleString(role);
			if (newRole != null) {
				roleString.append(delim + newRole);
				delim = " / ";
			}
		}
		// Store the role string.
		this.setStart(location, roleString.toString());
	}

	/**
	 * @return the representation of the role
	 *
	 * @param role	name of the role to represent
	 */
	protected abstract String roleString(String role);


	public static Map<String, ContigStarts> contigStartsMap(Genome genome, char strand,
			RoleMap roleMap) {
		Map<String, ContigStarts> retVal = new HashMap<String, ContigStarts>();
		// Create the contig entries.
		for (Contig contig : genome.getContigs()) {
			String contigID = contig.getId();
			ContigStarts startMap = (roleMap != null ? new Mapped(roleMap) : new UnMapped());
			retVal.put(contigID, startMap);
		}
		for (Feature feat : genome.getFeatures()) {
			Location floc = feat.getLocation();
			if (floc.getDir() == strand) {
				ContigStarts startMap = retVal.get(floc.getContigId());
				startMap.addFeature(feat);
			}
		}
		return retVal;
	}

	// SUBCLASSES

	/**
	 * Subclass in which role IDs are used.
	 */
	public static class Mapped extends ContigStarts {

		// FIELDS
		/** map for converting roles to role IDs */
		RoleMap roleMap;

		protected Mapped(RoleMap roles) {
			super();
			this.roleMap = roles;
		}

		@Override
		protected String roleString(String role) {
			Role roleObj = this.roleMap.getByName(role);
			String retVal = null;
			if (roleObj != null) retVal = roleObj.getId();
			return retVal;
		}

	}

	/**
	 * This class is used when there is no role map.  The roles are delimited by slashes.
	 */
	public static class UnMapped extends ContigStarts {

		protected UnMapped() {
			super();
		}

		@Override
		protected String roleString(String role) {
			return role;
		}

	}


}
