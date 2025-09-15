
package com.example.multidatasource.repository.c;

import com.example.multidatasource.model.c.EntityC;
import org.springframework.data.jpa.repository.JpaRepository;
import org.springframework.stereotype.Repository;

@Repository
public interface EntityCRepository extends JpaRepository<EntityC, Long> {
}
